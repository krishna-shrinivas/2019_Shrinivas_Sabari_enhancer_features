
# coding: utf-8

# In[1]:


import hoomd
import hoomd.md
import numpy
import render_figures
import fresnel
import gsd
import gsd.fl
import gsd.hoomd
import PIL
import os
import io
import math
import sys
import itertools
import csv
import matplotlib.pyplot as plt
import freud
import argparse
import random

device = fresnel.Device(mode='cpu');
path_tracer = fresnel.tracer.Path(device, 1920,1080)

blue = fresnel.color.linear([0.25,0.5,1])*0.9;
orange = fresnel.color.linear([1.0,0.714,0.169])*0.9

def render_sphere_frame(frame, height=None):

    if height is None:
        if hasattr(frame, 'configuration'):
            Ly = frame.configuration.box[1]
            height = Ly * math.sqrt(3)
        else:
            Ly = frame.box.Ly;
            height = Ly * math.sqrt(3)

    scene = fresnel.Scene(device)
    scene.lights = fresnel.light.cloudy();
    g = fresnel.geometry.Sphere(scene, position=frame.particles.position, radius=frame.particles.diameter*0.5)
    g.material = fresnel.material.Material(solid=1, color=blue, primitive_color_mix=1.0, specular=1.0, roughness=0.2)
    g.outline_width = 0.07
    scene.camera = fresnel.camera.orthographic(position=(height, height, height), look_at=(0,0,0), up=(0.2,1,0), height=height)

    g.color[frame.particles.typeid == 0] = blue;
    g.color[frame.particles.typeid == 1] = orange;
    g.color[frame.particles.typeid == 2] = [1,0,0];
    g.color[frame.particles.typeid == 3] = orange;

    scene.background_color = (1,1,1)

    return path_tracer.sample(scene, samples=64, light_samples=20)

def save_movie(frame_gen, gsd_file):
    f = gsd.fl.GSDFile(gsd_file, 'rb')
    t = gsd.hoomd.HOOMDTrajectory(f)

    a = frame_gen(t[0]);

    if tuple(map(int, (PIL.__version__.split(".")))) < (3,4,0):
        print("Warning! Movie display output requires pillow 3.4.0 or newer.")
        print("Older versions of pillow may only display the first frame.")

    im0 = PIL.Image.fromarray(a[:,:, 0:3], mode='RGB').convert("P", palette=PIL.Image.ADAPTIVE);
    ims = [];
    for f in t[1:]:
        a = frame_gen(f);
        im = PIL.Image.fromarray(a[:,:, 0:3], mode='RGB')
        im_p = im.quantize(palette=im0);
        ims.append(im_p)

    file_full_name = os.path.splitext(gsd_file)[0];
    path_name,file_name = os.path.split(file_full_name);
    f =path_name+'/Movies/'+file_name +'.gif';
    if not os.path.exists(os.path.dirname(f)):
        os.makedirs(os.path.dirname(f),exist_ok=True)

    im0.save(f, 'gif', save_all=True, append_images=ims, duration=1500, loop=0)

    return (f)



def random_initialize(position,L,N_min,N_end,min_gap):
    mol =N_min;
    while(mol<N_end):
        p = min_gap*numpy.floor(numpy.random.randint(low=  numpy.ceil(-L/2) ,high= numpy.floor(L/2), size=(1,3))/min_gap);
        if not any(numpy.equal(position,p).all(1)):
            position[mol,:]=p ;
            mol = mol+1;

    return position;

def generate_cubic_position(r):
    a = numpy.zeros((26,3));
    a[0:3,:]=-1*((numpy.array(list(set(itertools.permutations([0,0,r],3))))));
    a[3:6,:]= numpy.array((list(set(itertools.permutations([0,0,r],3)))));
    a[6:12,:]=numpy.array((list(set(itertools.permutations([0,r,-r],3)))));
    a[12:15,:]=numpy.array((list(set(itertools.permutations([0,-r,-r],3)))));
    a[15:18,:] = numpy.array((list(set(itertools.permutations([0,r,r],3)))));
    a[18:21,:] = numpy.array((list(set(itertools.permutations([-r,r,r],3)))));
    a[21:24,:] = numpy.array((list(set(itertools.permutations([-r,-r,r],3)))));
    a[24,:] = numpy.array((list(set(itertools.permutations([-r,-r,-r],3)))));
    a[25,:] = numpy.array((list(set(itertools.permutations([r,r,r],3)))));
    return a;

def generate_cubic_face_position(l,d):
    a = numpy.zeros((6,3));
    a[[0,2,4],:]=-1* numpy.array((list(set(itertools.permutations([0,0,0.5*(l-d)],3)))));
    a[[1,3,5],:]= numpy.array((list(set(itertools.permutations([0,0,0.5*(l-d)],3)))));
    return a;

class rdf_analyze:
    def __init__(self, system):
        self.system = system;
        self.rdf = freud.density.RDF(rmax=5.0, dr=0.01);

    def __call__(self, step):
        snap = system.take_snapshot();
        pos = snap.particles.position;
        N = system.N_total;
        box = freud.box.Box(snap.box.Lx, snap.box.Ly,snap.box.Lz);
        self.rdf.accumulate(box, pos[system.N_A:N,:],pos[system.N_A:N,:]);

class cluster_size_analyze:
    def __init__(self,system,cluster_size=None,cluster_rg=None):
        self.system = system;


    def __call__(self,step):
        snap = self.system.take_snapshot();
        N = self.system.N_total;
        pos =snap.particles.position[0:N,:];
        box = freud.box.Box(snap.box.Lx, snap.box.Ly,snap.box.Lz);
        cluster_g = freud.cluster.Cluster(box,rcut=1.4);
        cluster_g.computeClusters(pos)
        cluster_idx = cluster_g.getClusterIdx();
        cluster_prop = freud.cluster.ClusterProperties(box);
        cluster_prop.computeProperties(pos, cluster_idx)
        a = cluster_prop.getClusterSizes();
#        a_index = numpy.where(a==max(a))
##        rg = cluster_prop.getClusterG();
#        MI_rel = rg[a_index,:,:];
#        eig, eig_val = numpy.linalg.eig(MI_rel);

        return (max(a))

class cluster_rg_analyze:
    def __init__(self,system,cluster_size=None,cluster_rg=None):
        self.system = system;

    def __call__(self, step):
        snap = self.system.take_snapshot();
        N = self.system.N_total;
        pos = snap.particles.position[0:N,:];
        box = freud.box.Box(snap.box.Lx, snap.box.Ly,snap.box.Lz);
        cluster_g = freud.cluster.Cluster(box,rcut=1.4);
        cluster_g.computeClusters(pos)
        cluster_idx = cluster_g.getClusterIdx();
        cluster_prop = freud.cluster.ClusterProperties(box);
        cluster_prop.computeProperties(pos, cluster_idx)
        a = cluster_prop.getClusterSizes();
        a_index = numpy.where(a==max(a))
        rg = cluster_prop.getClusterG();
        MI_rel = rg[a_index,:,:];
#         print(MI_rel)
        eig, eig_val = numpy.linalg.eig(MI_rel);
#         print(numpy.sum(eig))
        return (numpy.sqrt(numpy.sum(eig))+0.5*self.system.diameter_list['B'])

#In[1.5]

def input_parse(filename):
    input_parameters  ={};
    with open(filename, 'r') as f:
        for line in f:
            line=line.strip();
            var_name,var_value = line.split(',');
            input_parameters[var_name] = float(var_value);
#            print("Value of {} is {}".format(var_name,input_parameters[var_name]))
    return input_parameters;

def params_parse(filename):
    param_list  =[];
    with open(filename, 'r') as f:
        param_name = f.readline().strip();
        for line in f:
            param_list.append(float(line.strip()))
    return (param_name,param_list);

# In[2]:
def generate_plots(file_output,f_cluster,file_PE_T_KE):
    from matplotlib import pyplot
    #get_ipython().magic('matplotlib inline')
    pyplot.figure(figsize=(4,2.2), dpi=140);
    data = numpy.genfromtxt(fname=f_cluster, skip_header=True);
    pyplot.plot(data[:,0]/1e5, data[:,1]);
    pyplot.xlabel('time step (* 10^5)');
    pyplot.ylabel('S(t)');
    pyplot.savefig(file_output+'_S.png',bbox_inches='tight')


    pyplot.figure(figsize=(4,2.2), dpi=140);
    pyplot.plot(data[:,0]/1e5, numpy.sqrt(data[:,2]) );
    pyplot.xlabel('time step (* 10^5)');
    pyplot.ylabel('Rg(t)');
    pyplot.savefig(file_output+'_Rg.png',bbox_inches='tight')



    pyplot.figure(figsize=(4,2.2), dpi=140);
    vol = 4/3 * numpy.pi * pow((data[:,2]),3)
    pyplot.plot(data[:,0]/1e5, numpy.divide(data[:,1],vol) );
    pyplot.xlabel('time step (* 10^5)');
    pyplot.ylabel('Rg(t)');
    pyplot.savefig(file_output+'_frac.png',bbox_inches='tight')


    # Use matplotlib to plot the potential energy vs time step.

    # In[15]:


    # #get_ipython().magic('matplotlib inline')
    # pyplot.figure(figsize=(4,2.2), dpi=140);
    # data = numpy.genfromtxt(fname=file_g_r_output, skip_header=False);
    # pyplot.plot(data[:,0], data[:,1]);
    # pyplot.xlabel('time step');
    # pyplot.ylabel('g(r)');


    # In[16]:


    data = numpy.genfromtxt(fname=file_PE_T_KE, skip_header=True);
    pyplot.figure(figsize=(4,2.2), dpi=140);
    pyplot.plot(data[:,0]/1e5, data[:,1]);
    pyplot.xlabel('time step (* 10^5)');
    pyplot.ylabel('potential_energy');
    pyplot.savefig(file_output+'_PE.png',bbox_inches='tight')


    # In[17]:


    pyplot.figure(figsize=(4,2.2), dpi=140);
    pyplot.plot(data[:,0]/1e5, data[:,2]);
    pyplot.xlabel('time step (* 10^5)');
    pyplot.ylabel('temperature');
    pyplot.savefig(file_output+'_T.png',bbox_inches='tight')
    # In[17]:


    pyplot.figure(figsize=(4,2.2), dpi=140);
    pyplot.plot(data[:,0]/1e5, data[:,3]);
    pyplot.xlabel('time step (* 10^5)');
    pyplot.ylabel('Kinetic energy');
    pyplot.savefig(file_output+'_KE.png',bbox_inches='tight')

    # In[18]:

    #
    # render_figures.render_sphere_frame(system.take_snapshot())

def write_input_params(file_output,input_params):
    with open(file_output,'w') as f:
        for key in input_params:
            f.write(key+','+str(input_params[key])+'\n');
    f.close()


def run_MD_simulation(input_file,output_file,params_file):

    input_params = input_parse(input_file);

    if params_file is not None:
        (param_name,param_list) = params_parse(params_file)
    else:
        param_name = 'N_A';
        param_list = [];
        param_list.append(int(input_params[param_name]));

    for param in param_list:
        hoomd.context.initialize("");

        input_params[param_name] = param;

        #Particle counts of species A,B, and C
        N_A = int(input_params['N_A']);
        N_B =int(input_params['N_B']);
        N_C = int(input_params['N_C']);


        # Interaction energy between the various types of particles
        base_interaction = 0.05;

        weak_interaction_BC = input_params['weak_interaction_BC'];
        weak_interaction_BB = input_params['weak_interaction_BB'];
        weak_interaction_CC = input_params['weak_interaction_CC'];


        strong_interaction_AB =input_params['strong_interaction_AB'];
        strong_interaction_AC =input_params['strong_interaction_AC'];


        seed_position = int(input_params['seed_position']);
        # Size of cubic box, and temperature of system
        L =int(input_params['L']);
        L_A_box = int(input_params['L_A_box']);
        T = input_params['T'];

        # Diameter of various particles, and patches
        particle_types = ['A', 'B','C','ACS','ABS','BDNA','CDNA','AW','GS'];
        diameter_list = {}
        for ptype in particle_types:
            diameter_list[ptype] = input_params['d_'+ptype]

        N_bs_AB = int(input_params['N_bs_AB']);
        N_bs_AC = int(input_params['N_bs_AC']);

        log_step = int(input_params['log_step']);
        log_trajectory_step = int(input_params['log_trajectory_step']);
        t_end =int(input_params['t_end']);
        seed = int(input_params['seed']);
        seed = random.randint(1,1000000000)

        drug_test = int(input_params['drug_test']);
        t_post_drug = int(input_params['t_post_drug']);
        p_drug = float(input_params['p_drug'])



        meta_data = str(N_A)+'_'+str(N_B)+ '_'+str(N_C)+'_'+str(L)+'_'+str(seed_position);

        file_output = 'Output/'+meta_data+'/Data/Observables/'+param_name+'_'+str(param)+'/'+output_file+'_'+ str(seed);
        trajectory_file = 'Output/'+meta_data+'/Data/Trajectory/'+param_name+'_'+str(param)+'/'+output_file+'_'+ str(seed)+'.gsd';
        figure_output =  'Output/'+meta_data+'/Data/Figure/'+param_name+'_'+str(param)+'/'+output_file+'_'+ str(seed);

        if not os.path.exists(os.path.dirname(file_output)):
            os.makedirs(os.path.dirname(file_output),exist_ok=True)

        if not os.path.exists(os.path.dirname(trajectory_file)):
            os.makedirs(os.path.dirname(trajectory_file),exist_ok=True)

        if not os.path.exists(os.path.dirname(figure_output)):
            os.makedirs(os.path.dirname(figure_output),exist_ok=True)




        file_input_params = file_output + '_input_params.text';
        write_input_params(file_input_params,input_params);
        N_total =N_A+N_B+N_C;

        snapshot = hoomd.data.make_snapshot(N=N_total, box=hoomd.data.boxdim(L=L),  particle_types=particle_types);

        #Assign particle type ID's
        snapshot.particles.typeid[0:N_A] = 0;
        snapshot.particles.typeid[N_A:N_A+N_B] = 1;
        snapshot.particles.typeid[N_B+N_A:N_total] = 2;

        #Assign particle position
        position = numpy.zeros((N_total,3));
        numpy.random.seed(seed=seed_position);
        position =random_initialize(position,L_A_box,0,N_A,1)
        position =random_initialize(position,L,N_A,N_total,1)
        # print(position)
        snapshot.particles.position[:] = position;



        # Set particle moment of inertia, diameter, and velocity
        snapshot.particles.moment_inertia[N_A:N_total,:] = numpy.ones((N_total-N_A,3))
        snapshot.particles.mass[0:N_total] = 1;
        snapshot.particles.diameter[0:N_A] = diameter_list['A'];
        snapshot.particles.diameter[N_A:N_A+N_B] = diameter_list['B'];
        snapshot.particles.diameter[N_A+N_B:] = diameter_list['C'];
        snapshot.particles.velocity[N_A:] = numpy.random.normal(0.0,
                                        numpy.sqrt(T / 1.0), [snapshot.particles.N-N_A, 3]);




        #Intialize system with snapshot and store some system properties
        system=hoomd.init.read_snapshot(snapshot)
        system.N_A = N_A;
        system.N_total = N_total;
        system.N_bs_A = N_bs_AB;
        system.N_bs_AC = N_bs_AC;

        system.diameter_list = diameter_list;
        print(system.particles.types)

        render_figures.render_sphere_frame(system.take_snapshot())




        # Set up the rigid structure of the bodies
        rigid = hoomd.md.constrain.rigid();

        # Generates a function that constructs a cubic particle
        pos = generate_cubic_position(diameter_list['A']);

        #Generates positions of mid faces, scaled by binding site size
        face_pos = (generate_cubic_face_position(3*diameter_list['A'],diameter_list['ABS']))


        pos = numpy.vstack([pos,face_pos])

        diam = numpy.vstack([diameter_list['AW']*numpy.ones((26,1)),diameter_list['ABS']*numpy.ones((6,1))]);
        N_bs = system.N_bs_A;
        N_bs_AC = system.N_bs_AC;

        rigid.set_param('A',
                        types=['AW']*26+['ABS']*N_bs+['ACS']*N_bs_AC,
                        positions=pos[0:N_bs_AC+N_bs+26,:],orientations=None, charges=None, diameters=diam[0:N_bs_AC+N_bs+26]);


        rigid.set_param('B',
                        types=['BDNA']*1,
                        positions=[(0,0,-0.5*(diameter_list['B'])+0.5*diameter_list['BDNA'])],orientations=None, charges=None, diameters=[diameter_list['BDNA']]);

        rigid.set_param('C',
                        types=['CDNA']*1,
                        positions=[(0,0,-0.5*(diameter_list['C'])+0.5*diameter_list['CDNA'])],orientations=None, charges=None, diameters=[diameter_list['CDNA']]);



        rigid.create_bodies()
        render_figures.render_sphere_frame(system.take_snapshot())





        # Estimate the actual volume fraction

        vol_A = N_A*pow(3*diameter_list['A'],3)/float(pow(L,3))
        vol_B = N_B*4/3*numpy.pi*pow(diameter_list['B']/2,3)/float(pow(L,3));
        vol_C = N_C*4/3*numpy.pi*pow(diameter_list['C']/2,3)/float(pow(L,3));

        vol_fraction = vol_A+vol_B+vol_C;
        print("The volume fraction is {}".format(vol_fraction))
        print("The individual volume fraction are B = {} and C= {}".format(vol_B,vol_C))




        hoomd.md.integrate.mode_standard(dt=0.005,aniso=True);
        groupB = hoomd.group.type(name='groupB', type='B')
        groupC = hoomd.group.type(name='groupC', type='C')
        groupBC = hoomd.group.union(name="bc-particles", a=groupB, b=groupC)
        integrator= hoomd.md.integrate.langevin(group=groupBC, kT=T, seed=seed);




        all = hoomd.group.all()
        file_PE_T_KE = file_output +'_PE_T_KE.log'
        print((file_PE_T_KE))
        hoomd.analyze.log(filename=file_PE_T_KE,
                          quantities=['potential_energy', 'temperature','kinetic_energy'],
                          period=log_step,
                          overwrite=True);
        hoomd.dump.gsd(trajectory_file, period=log_trajectory_step, group=all, overwrite=True);




        nl = hoomd.md.nlist.cell()



        def lj_shifted(r, rmin, rmax, epsilon,d1,d2):
            sigma = 0.5*(d1+d2);
            r_max = 2.5*sigma;
            V_shift = 4 * epsilon * ( (sigma / r_max)**12 - (sigma / r_max)**6);
            V = 4 * epsilon * ( (sigma / r)**12 - (sigma / r)**6) - V_shift;
            F = 4 * epsilon / r * ( 12 * (sigma / r)**12 - 6 * (sigma / r)**6);
            return (V, F)

        table = hoomd.md.pair.table(width=1000,nlist=nl);

        ci = 0;

        for typei in system.particles.types:
            cj=0;
            for typej in system.particles.types:
                table.pair_coeff.set(typei, typej, func=lj_shifted, rmin=0.1, rmax=1.25*(system.diameter_list[typei]+system.diameter_list[typej]), coeff=dict(epsilon=base_interaction, d1=system.diameter_list[typei],d2=system.diameter_list[typej]));
                cj  = cj+1;
            ci = ci+1;

        table.pair_coeff.set('ABS', 'BDNA', func=lj_shifted, rmin=0.1, rmax=1.25*(system.diameter_list['ABS']+system.diameter_list['BDNA']), coeff=dict(epsilon=strong_interaction_AB, d1=system.diameter_list['ABS'],d2=system.diameter_list['BDNA']))
        table.pair_coeff.set('ACS', 'CDNA', func=lj_shifted, rmin=0.1, rmax=1.25*(system.diameter_list['ACS']+system.diameter_list['CDNA']), coeff=dict(epsilon=strong_interaction_AC, d1=system.diameter_list['ACS'],d2=system.diameter_list['CDNA']))

        table.pair_coeff.set('C', 'C', func=lj_shifted, rmin=0.1, rmax=1.25*(system.diameter_list['C']+system.diameter_list['C']), coeff=dict(epsilon=weak_interaction_CC, d1=system.diameter_list['C'],d2=system.diameter_list['C']))
        table.pair_coeff.set('B', 'C', func=lj_shifted, rmin=0.1, rmax=1.25*(system.diameter_list['B']+system.diameter_list['C']), coeff=dict(epsilon=weak_interaction_BC, d1=system.diameter_list['B'],d2=system.diameter_list['C']))
        table.pair_coeff.set('B', 'B', func=lj_shifted, rmin=0.1, rmax=1.25*(system.diameter_list['B']+system.diameter_list['B']), coeff=dict(epsilon=weak_interaction_BB, d1=system.diameter_list['B'],d2=system.diameter_list['B']))



    #    analyzer = rdf_analyze(system);
    #    hoomd.analyze.callback(analyzer, period=log_step, phase= 0);
    #


    #    In[10]:


        f_cluster = file_output+'_cluster_size.dat';
        logger = hoomd.analyze.log(filename=f_cluster, quantities=['cluster_size','cluster_rg'], period=log_step,overwrite=True)
        a = cluster_size_analyze(system)
        logger.register_callback('cluster_size', cluster_size_analyze(system))
        logger.register_callback('cluster_rg', cluster_rg_analyze(system))


        # Take 10,000 steps forward in time.12

     #In[ ]:


        hoomd.run(t_end);


        #In[12]:

        if drug_test:
            N_types = len(system.particles.types);
            print(N_types)
            groupCDNA = hoomd.group.type(name="cdna-particles", type='CDNA')
            for p in groupCDNA:
                if p_drug >= numpy.random.random_sample():
                    p.type = 'GS'

            hoomd.run(t_post_drug);


        #In[13]:
        generate_plots(figure_output,f_cluster,file_PE_T_KE)
        save_movie(render_sphere_frame, trajectory_file);
    #    file_g_r_output = file_output+'_rdf.log'
    #    with open(file_g_r_output, 'w') as f:
    #        writer = csv.writer(f, delimiter='\t')
    #        writer.writerows(zip(analyzer.rdf.getR(),analyzer.rdf.getRDF()))

        print("Simulation is complete now!");
# ## Examine the output

# In[11]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Take input and output filename and run MD simulations')
    parser.add_argument('--i',help="Name of the input_file", required = True);
    parser.add_argument('--o',help="Name of output file", required = True);
    parser.add_argument('--p',help="List of changing parameters", required = False)
    args = parser.parse_args();
    # params_file = args.param_list;
    # if params_file is not None:
    #     (param_name,param_list) = params_parse(params_file)

    run_MD_simulation(args.i,args.o,args.p);

# In[14]:


# Examine how the system configuration evolves over time. [render_figures](render_figures.py) is a helper script that builds animated gifs from trajectory files and system snapshots. It is part of the [hoomd-examples](https://bitbucket.org/glotzer/hoomd-examples) repository and designed only to render these examples.

# In[19]:


#import render_figures
#render_figures.display_movie(render_figures.render_sphere_frame, trajectory_file);
#
