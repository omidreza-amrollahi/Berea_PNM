import numpy as np
import pandas as pd
import porespy as ps
from porespy.networks import regions_to_network, add_boundary_regions
from porespy.networks import _net_dict
from porespy.networks import label_boundary_cells
from porespy.tools import pad_faces
from porespy.tools import make_contiguous
from porespy.metrics import region_surface_areas, region_interface_areas
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import openpnm as op
import scipy as sp
from skimage import io
from skimage.morphology import binary_dilation
import streamlit as st

from openpnm.network import GenericNetwork
from openpnm.geometry import Imported
from openpnm.io import GenericIO
from openpnm.models import physics as mods


def create_workspace():
    ws = op.Workspace()
    ws.clear()
    return ws

def get_porosity(file):
    im = io.imread(file)
    im = im > 1
    imtype=im.view()
    im = sp.array(im, dtype=bool)
    im = ~im 
    im = im.T
    st.write(f"Initial Porosity: {ps.metrics.porosity(im).round(3)}")
    return im

def plot_porosity_profile(im):
    ps.visualization.set_mpl_style()
    voxel_size = 1.843803  # microns/voxel
    x_profile = ps.metrics.porosity_profile(im, 0)
    y_profile = ps.metrics.porosity_profile(im, 1)
    z_profile = ps.metrics.porosity_profile(im, 2)
    fig, ax = plt.subplots()
    ax.plot(np.linspace(0, im.shape[0]*voxel_size, im.shape[0]), x_profile, 'b-', label='yz-plane', alpha=0.5)
    ax.plot(np.linspace(0, im.shape[1]*voxel_size, im.shape[1]), y_profile, 'r-', label='xz-plane', alpha=0.5)
    ax.plot(np.linspace(0, im.shape[2]*voxel_size, im.shape[2]), z_profile, 'g-', label='xy-plane', alpha=0.5)
    ax.set_ylabel('Porosity of slice')
    ax.set_xlabel('Position of slice along given axis')
    ax.legend()
    st.pyplot(fig)

@st.cache()
def get_inscribed_sphere(regions):
    regions_temp = regions.regions*regions.im
    props = ps.metrics.regionprops_3D(regions_temp)
    sph = ps.metrics.props_to_image(regionprops=props, shape=im.shape, prop='inscribed_sphere')
    return sph  
    
def plot_inscribed_sphere(sph, slice=1):
    fig, ax = plt.subplots()
    ax.imshow(sph[slice,:,:] + 0.5*(~im[slice,:,:]) , cmap=plt.cm.inferno)
    st.pyplot(fig)

@st.cache()
def snow_partitioning(im, r_max, sigma, return_all=True):
    regions = ps.filters.snow_partitioning(im=im,  r_max=r_max, sigma=sigma, return_all=return_all)
    return regions

def plot_snow_partitioning(regions, slice=1):
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=[10, 10])
    ax1.imshow(regions.im[slice,:,:], origin='lower')
    ax2.imshow(regions.dt[slice,:,:], origin='lower')
    dt_peak = regions.dt.copy()
    peaks_dilated = binary_dilation(regions.peaks > 0)
    dt_peak[peaks_dilated > 0] = np.nan
    cmap = cm.viridis
    cmap.set_bad('red', 1.)
    ax3.imshow(dt_peak[slice,:,:], origin='lower', cmap=cmap)
    ax4.imshow(regions.regions[slice,:,:], origin='lower')
    st.pyplot(fig)
    st.write(f"Number of regions: {np.unique(regions.regions).size}")

@st.cache(allow_output_mutation=True)
def extract_network(regions, boundary_faces, voxel_size, marching_cubes_area=False):
    im = regions.im
    dt = regions.dt
    regions = regions.regions
    b_num = sp.amax(regions)
    # -------------------------------------------------------------------------
    # Boundary Conditions
    regions = add_boundary_regions(regions=regions, faces=boundary_faces)
    # -------------------------------------------------------------------------
    # Padding distance transform and image to extract geometrical properties
    dt = pad_faces(im=dt, faces=boundary_faces)
    im = pad_faces(im=im, faces=boundary_faces)
    regions = regions*im
    regions = make_contiguous(regions)
    # -------------------------------------------------------------------------
    # Extract void and throat information from image
    net = regions_to_network(im=regions, dt=dt, voxel_size=voxel_size)
    # -------------------------------------------------------------------------
    # Extract marching cube surface area and interfacial area of regions
    if marching_cubes_area:
        areas = region_surface_areas(regions=regions)
        interface_area = region_interface_areas(regions=regions, areas=areas,
                                                voxel_size=voxel_size)
        net['pore.surface_area'] = areas * voxel_size**2
        net['throat.area'] = interface_area.area
    # -------------------------------------------------------------------------
    # Find void to void connections of boundary and internal voids
    boundary_labels = net['pore.label'] > b_num
    loc1 = net['throat.conns'][:, 0] < b_num
    loc2 = net['throat.conns'][:, 1] >= b_num
    pore_labels = net['pore.label'] <= b_num
    loc3 = net['throat.conns'][:, 0] < b_num
    loc4 = net['throat.conns'][:, 1] < b_num
    net['pore.boundary'] = boundary_labels
    net['throat.boundary'] = loc1 * loc2
    net['pore.internal'] = pore_labels
    net['throat.internal'] = loc3 * loc4
    # -------------------------------------------------------------------------
    # label boundary cells
    net = label_boundary_cells(network=net, boundary_faces=boundary_faces)
    # -------------------------------------------------------------------------
    # assign out values to dummy dict

    temp = _net_dict(net)
    temp.im = im.copy()
    temp.dt = dt
    temp.regions = regions
    net = temp
    return net

def get_pore_network(net):
    network = GenericNetwork(project=None)
    network = GenericIO._update_network(network=network, net=net)
    Imported(network=network, settings={})
    prj = network.project
    pn = prj.network
    geom = prj.geometries()['geo_01']
    net_health = pn.check_network_health()
    op.topotools.trim(network=pn, pores=net_health["trim_pores"])
    return pn, prj, geom

def plot_network(pn, face):
    fig, ax = plt.subplots()
    op.topotools.plot_connections(network=pn, alpha=0.8, color='grey', ax=ax)
    op.topotools.plot_coordinates(network=pn, ax=ax, color='b', markersize=10)
    op.topotools.plot_coordinates(network=pn, pores=pn.pores(face), ax=ax, color='r', markersize=100)
    st.pyplot(fig)
    
def modify_network(pn, geom):
    del geom['pore.area']
    del geom['pore.volume']
    del geom['throat.conduit_lengths.pore1']
    del geom['throat.conduit_lengths.pore2']
    del geom['throat.conduit_lengths.throat']
    del geom['throat.endpoints.tail']
    del geom['throat.endpoints.head']


    geom['throat.volume'] = np.zeros(geom['throat.volume'].shape)
    geom.add_model(propname='pore.volume',
                    model=op.models.geometry.pore_volume.sphere)
    geom['pore.volume'][pn.pores('boundary')] = 0
    geom.add_model(propname='throat.endpoints',
                    model=op.models.geometry.throat_endpoints.spherical_pores)
    geom.add_model(propname='pore.area',
                    model=op.models.geometry.pore_area.sphere)
    geom.add_model(propname='throat.conduit_lengths',
                    model=op.models.geometry.throat_length.conduit_lengths)
    return geom

def get_phase(pn, name, throatCA, poreCA, throatSigma, poreSigma, throatmu, poremu):
    phase = op.phases.GenericPhase(network=pn, name=name)
    phase['throat.contact_angle'] = throatCA
    phase['pore.contact_angle'] = poreCA
    phase['throat.surface_tension'] = throatSigma
    phase['pore.surface_tension'] = poreSigma
    phase['throat.viscosity'] = throatmu
    phase['pore.viscosity'] = poremu
    return phase

def add_physics(pn, phase, geom):
    phys = op.physics.GenericPhysics(network=pn, phase=phase, geometry=geom)
    phys.add_model(propname='throat.entry_pressure',
                    model=mods.capillary_pressure.washburn)
    phys['pore.entry_pressure'] = 0
    phys.add_model(propname='throat.flow_shape_factors',
               model=mods.flow_shape_factors.conical_frustum_and_stick)
    phys.add_model(propname='throat.hydraulic_conductance',
                model=mods.hydraulic_conductance.hagen_poiseuille)
    return phys

def plot_ip(alg_ip):
    data_ip = alg_ip.get_intrusion_data()
    fig, ax = plt.subplots()
    ax.semilogy(data_ip.S_tot, data_ip.Pcap)
    ax.set_xlabel('Invading Phase Saturation')
    ax.set_ylabel('Capillary Pressure (bars)')
    ax.grid(True)
    st.pyplot(fig)

@st.cache()
def perform_ip(pn, oil, face):
    alg_ip = op.algorithms.InvasionPercolation(network=pn)
    alg_ip.setup(phase=oil)
    alg_ip.set_inlets(pores=pn.pores(face))
    alg_ip.run()
    # alg_ip.apply_trapping(outlets=pn.pores(['opposite of face']))
    oil.update(alg_ip.results())
    data = alg_ip.get_intrusion_data()
    return alg_ip, pd.DataFrame(data, index=["capillary pressure (bar)", "invading phase saturation"]).transpose().to_csv().encode('utf-8')

def calculate_permeability(pn, oil):
    Pin = 1.0e5
    Pout = 0.0


    K_single_phase = [None,None,None]
    bounds = [ ['top', 'bottom'], ['left', 'right'],['front', 'back']]
    [amax, bmax, cmax] = np.max(pn['pore.coords'], axis=0)
    [amin, bmin, cmin] = np.min(pn['pore.coords'], axis=0)
    lx = amax-amin
    ly = bmax-bmin
    lz = cmax-cmin
    da = lx*ly
    dl = lz


    def top_b(lx,ly,lz):
        da = lx*ly
        dl = lz
        res_2=[da,dl]
        return res_2

    def left_r(lx,ly,lz):
        da = lx*lz
        dl = ly
        res_2=[da,dl]
        return res_2

    def front_b(lx,ly,lz):
        da = ly*lz
        dl = lx
        res_2=[da,dl]
        return res_2

    options = {0 : top_b(lx,ly,lz),1 : left_r(lx,ly,lz),2 : front_b(lx,ly,lz)}
    for bound_increment in range(len(bounds)):
        BC1_pores = pn.pores(labels=bounds[bound_increment][0])
        BC2_pores = pn.pores(labels=bounds[bound_increment][1])
        [da,dl]=options[bound_increment]
        
        # Permeability - water
        sf = op.algorithms.StokesFlow(network=pn, phase=oil)
        sf.setup(conductance='throat.hydraulic_conductance')
        sf.set_value_BC(pores=BC1_pores, values=Pin)
        sf.set_value_BC(pores=BC2_pores, values=Pout)
        sf.run()
        K_single_phase[bound_increment] = sf.calc_effective_permeability(domain_area=da,
                                                                            domain_length=dl,
                                                                            inlets=BC1_pores,
                                                                            outlets=BC2_pores)
        fig, ax = plt.subplots()
        psm = op.topotools.plot_coordinates(pn, color=sf['pore.pressure'], 
                                        size_by=pn['pore.diameter'], 
                                        markersize=50,
                                        ax = ax)
        fig.colorbar(psm, ax=ax, location='left')
        st.pyplot(fig)
    
    st.write(f"The permeability coefficient from top to bottom is: {K_single_phase[0][0].round(3)}")
    st.write(f"The permeability coefficient from left to right is: {K_single_phase[1][0].round(3)}")
    st.write(f"The permeability coefficient from front to back is: {K_single_phase[2][0].round(3)}")

@st.cache()
def calculate_relative_permeability(pn, invading_phase_name, defending_phase_name):
    rp = op.algorithms.metrics.RelativePermeability(network=pn)
    rp.setup(invading_phase=invading_phase_name, defending_phase=defending_phase_name,
            invasion_sequence='invasion_sequence')
    rp.run()
    results=rp.get_Kr_data()
    return pd.DataFrame(results).to_csv().encode('utf-8'), rp

st.title("Pore network modelling of a porous medium using OpenPNM")
st.write("This app is a demonstration of the use of OpenPNM for modelling a porous medium. The app allows you to upload a tiff image of a porous medium and then perform a series of calculations on the image")
st.write("Aknowledgegement: This app was developed by [Omidreza Amrollahinasab](https://github.com/omidreza-amrollahi) from [Department of Petroleum Engineering Leoben](dpe.at)")

file = st.file_uploader("Upload a tiff image", type=["tiff", "tif"])
if file is not None:
    ws = create_workspace()
    im = get_porosity(file)
    
    if st.button("Plot Porosity Profile"):
        plot_porosity_profile(im)
    
    voxel_size = st.number_input("Enter voxel size in microns", value=1.843803)
    r_max = st.number_input("Enter maximum pore size in microns", value=5)
    sigma = st.number_input("Enter sigma value", value=0.35)
    boundary_faces=['top', 'bottom', 'left', 'right', 'front', 'back']
    marching_cubes_area=False
    
    regions = snow_partitioning(im, r_max, sigma, return_all=True)
    slice = st.number_input("Enter slice number to plot along x axis", value=1)
    plot_snow_partitioning(regions, slice=slice)
    
    sph = get_inscribed_sphere(regions)
    plot_inscribed_sphere(sph, slice=slice)
    
    net = extract_network(regions, boundary_faces, voxel_size, marching_cubes_area)
    pn, prj, geom = get_pore_network(net)

    face = st.selectbox("Select a face", boundary_faces)
    plot_network(pn, face)

    geom = modify_network(pn, geom)
    
    p1_throatCA = st.number_input("Enter invading phase throat contact angle in degrees", value=180)
    p1_poreCA = st.number_input("Enter invading phase pore contact angle in degrees", value=180)
    p1_throatSigma = st.number_input("Enter invading phase throat surface tension in N/m", value=20e-3)
    p1_poreSigma = st.number_input("Enter invading phase pore surface tension in N/m", value=20e-3)
    p1_throatmu = st.number_input("Enter invading phase throat viscosity in Pa.s", value=0.8e-3)
    p1_poremu = st.number_input("Enter invading phase pore viscosity in Pa.s", value=0.8e-3)
    oil = get_phase(pn, 'oil', p1_throatCA, p1_poreCA, p1_throatSigma, p1_poreSigma, p1_throatmu, p1_poremu)
    add_physics(pn, oil, geom)
    alg_ip, indtrusion_data = perform_ip(pn, oil, face)
    plot_ip(alg_ip)
    st.download_button(
        label="Download capillary pressure as CSV",
        data=indtrusion_data,
        file_name='pc.csv',
        mime='text/csv',
    )
    with st.spinner('Calculating Permeability...'):
        calculate_permeability(pn, oil)
    
    p2_throatCA = st.number_input("Enter defending phase throat contact angle in degrees", value=0)
    p2_poreCA = st.number_input("Enter defending phase pore contact angle in degrees", value=0)
    p2_throatSigma = st.number_input("Enter defending phase throat surface tension in N/m", value=20e-3)
    p2_poreSigma = st.number_input("Enter defending phase pore surface tension in N/m", value=20e-3)
    p2_throatmu = st.number_input("Enter defending phase throat viscosity in Pa.s", value=9.4e-3)
    p2_poremu = st.number_input("Enter defending phase pore viscosity in Pa.s", value=9.4e-3)
    water = get_phase(pn, 'water', p2_throatCA, p2_poreCA, p2_throatSigma, p2_poreSigma, p2_throatmu, p2_poremu)
    add_physics(pn, water, geom)
    
    kr, rp = calculate_relative_permeability(pn, 'oil', 'water')
    fig, ax = plt.subplots()
    ax.set_xlim([0, 1])
    rp.plot_Kr_curves(fig)
    st.pyplot(fig)
    
    st.download_button(
        label="Download relative permeability as CSV",
        data=kr,
        file_name='kr.csv',
        mime='text/csv',
    )
        