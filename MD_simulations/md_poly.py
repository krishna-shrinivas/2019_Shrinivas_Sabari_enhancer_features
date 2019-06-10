"""
Created on Tue Mar 13 18:57:58 2018

@author: krishna
"""


"""
We are importing all the libraries that are required for running utils and
functions in this code. The primary dependencies include:

    HOOMD-blue
    Numpy
    Random
    argparse
    OS
    csv
    fresnel
    freud
    itertools
"""
import hoomd
import hoomd.md
import numpy
import os
import argparse
import random


from graphics.plots import *
from utils.intialize import *
from utils.input_parse import *
from utils.observables import *
from utils.potentials import *



def run_MD_simulation(input_file,output_file,pt_file,params_file,energy_file,outfolder,dest_folder):

    """
    Inputs to the function are the following:
        input_file --> path to file with input parameters
            Input file is a list of all input parameters including the number
            of chains, sequence of different types, size of simulation of box,
            size of simulation box for chain A initialization, length of simulation,
            dt for simulation, bond strength, bjeurrum length, screening parameter,
            drug characterization, and p_drug


        output_file --> prefix to trajectory file
            The prefix to the output file will contain the name of the prefix
            will be suffixed by the unique PBS ID


        pt_file --> path to file with monomer information
            Monomer information file which contains information on unique
            monomer types, their mass, their charge, and their diameter. This
            should include details for the patches which will model the unique
            DNA patches, unique BS pairs.


        params_file --> path to file with list of parameter sweep regimes
            Parameter file contains the unique name of  the parameter to sweep
            or the pairwise interactionns (for e.g. E,C) on the first line,
            and every subsequent line contains a value of the parameter. The code
            does not check that the param list is non-repetitive, so it'll redo
            simulations across multiple simulation regimes. In the absence
            of a parameter sweep flag being passed, the program assumes that
            param variable takes a value None.


        energy_file --> path to file with list of pairwise interactions
            For all pairwise LJ parameters not defined, a default
            value of 0.01 is assumed for E.

        outfolder --> Prefix to output data files
        dest_folder --> Destination folder which Output/out_folder will write

    """

#    input_file = 'input_parameters.text'
    
    if dest_folder is None:
        dest_folder = '';
    input_params = input_parse(input_file);
    pt_info=read_particle_type_information(pt_file);
    interaction_list = energetics_parse(energy_file);
    #read input parameters

#    pt_file = 'input_particle_information.text'
#    energy_file = 'energetics_parameters.text'

#
#    params_file = None;
#    output_file = 'traj'
#
    if params_file is not None:
        (param_name,param_list) = params_parse(params_file)
    else:
        param_name = 'seed_position';
        param_list = [];
        param_list.append(int(input_params[param_name]));

    if outfolder is None:
        outfolder = 'Default';

    for param in param_list:


        # Initialize hoomd simulation run
        # for each parameter value

        hoomd.context.initialize("");

        #   This checks if the param_name contains
        #   a ',' which only happens if the param
        #   is a energy of interaction
        if param_name.find(',') > -1:
            (param_name,param_sub_name) = param_name.split(',');

        if param_name in input_params.keys():
            input_params[param_name] = param;
        elif param_name in pt_info.keys():
            pt_info[param_name][param_sub_name] = param;
        else:
            interaction_list[(param_name,param_sub_name)]=param;
            interaction_list[(param_sub_name,param_name)]=param;

        #Time-steps of simulation
        dt = float(input_params['dt']);

        #Particle counts of species A,B, and C
        N_A = int(input_params['N_A']);
        N_B =int(input_params['N_B']);
        N_C = int(input_params['N_C']);

        #Sequences of A,B, and C
        seq_A = (input_params['seq_A']);
        seq_B =(input_params['seq_B']);
        seq_C = (input_params['seq_C']);

        #Temperature of A
        T_A = float(input_params['T_A']);

        #Monomer information
        particle_types = pt_info['type'];
        diameter_list = pt_info['diameter'];
        charge_list = pt_info['charge'];
        is_active_DNA_list = pt_info['is_active_DNA'];

        # Interaction energy between the various types of particles
        base_interaction = 0.01;

        for typei in particle_types:
            for typej in particle_types:
                if (typei,typej) not in interaction_list.keys():
                    interaction_list[(typei,typej)] = base_interaction;
                    interaction_list[(typej,typei )] = base_interaction;



        #movie-Flag information
        movie_flag = int(input_params['movie_flag']);
        plot_flag = int(input_params['plot_flag']);

        seed_position = int(input_params['seed_position']);


        # Size of cubic box, and temperature of system
        L =int(input_params['L']);
        L_A_box = int(input_params['L_A_box']);
        T = float(input_params['T']);

        # Number of patches of different types on A
        N_bs_AB = int(input_params['N_bs_AB']);
        N_bs_AC = int(input_params['N_bs_AC']);

        # Bond_information
        k = float(input_params['K']);

        #Charged media information
        # inverse of screening length and bjeurrum length
        kappa = float(input_params['kappa']);
        lb = float(input_params['lb']);

        #Damping coeffiecient for A-type molecules
        gamma_A =float(input_params['gamma_A']);

        # Simulation settings -
        # frequency of logging observables, trajectories, length of simulation
        # seed for the langevin dynamics is randomly generated

        log_step = int(input_params['log_step']);
        log_trajectory_step = int(input_params['log_trajectory_step']);
        t_end =float(input_params['t_end']);
        seed = random.randint(1,1000000000)
        print(seed)

        # Parameters for setting the JQ1 interactions
        # True or False for perturbation, time to track after drug, and dose of drug
        drug_test = int(input_params['drug_test']);
        t_post_drug = float(input_params['t_post_drug']);
        p_drug = float(input_params['p_drug'])
        
        
        if 'warm_up_flag' in input_params.keys():
            warm_up_flag = int(input_params['warm_up_flag']);

        else:
            warm_up_flag =0;
            
        if 'warm_up_steps' in input_params.keys():
            warm_up_steps = int(input_params['warm_up_steps']);
        else:
            warm_up_steps = 20000;
            print("Values are {},{}".format(warm_up_flag,warm_up_steps))

        meta_data = str(N_A)+'_'+str(N_B)+ '_'+str(N_C)+'_'+str(L)+'_'+str(seed_position)+'_'+str(seq_A)+'_'+str(seq_B)+'_'+str(seq_C);

        file_output = dest_folder + 'Output/' +outfolder+'/' +meta_data+'/Data/Observables/'+param_name+'_'+str(param)+'/'+output_file+'_'+ str(seed);
        trajectory_file = dest_folder + 'Output/' + outfolder+'/'+meta_data+'/Data/Trajectory/'+param_name+'_'+str(param)+'/'+output_file+'_'+ str(seed)+'.gsd';
        figure_output =  dest_folder + 'Output/'+outfolder+'/' +meta_data+'/Data/Figure/'+param_name+'_'+str(param)+'/'+output_file+'_'+ str(seed);

        if not os.path.exists(os.path.dirname(file_output)):
            os.makedirs(os.path.dirname(file_output),exist_ok=True)

        if not os.path.exists(os.path.dirname(trajectory_file)):
            os.makedirs(os.path.dirname(trajectory_file),exist_ok=True)

        if not os.path.exists(os.path.dirname(figure_output)):
            os.makedirs(os.path.dirname(figure_output),exist_ok=True)




        file_input_params = file_output + '_input_params.text';
        write_input_params(file_input_params,input_params);
        header_sep = {};
        header_sep['%']='interaction';
        write_input_params(file_input_params,header_sep);
        write_input_params(file_input_params,interaction_list);
        header_sep['%']='diameter'
        write_input_params(file_input_params,header_sep);
        write_input_params(file_input_params,diameter_list);
        header_sep['%']='charge'
        write_input_params(file_input_params,header_sep);
        write_input_params(file_input_params,charge_list);
        header_sep['%']='mass'
        write_input_params(file_input_params,header_sep);
        write_input_params(file_input_params,pt_info['mass']);



        # Calculate the total number of molecules, the total number of monomers, and bonds
        N_total =N_A*len(seq_A)+N_B*len(seq_B)+N_C*len(seq_C);
        N_mol = N_A+N_B+N_C;
        N_bonds = N_A*(len(seq_A)-1)+N_B*(len(seq_B)-1)+N_C*(len(seq_C)-1);

        # Assign each monomer to it's 'parent' polymer
        # and set up the simulation box
        position = numpy.zeros((N_total,3));

        pos_params ={};


        particle_info = {};
        particle_info['typeid'] = numpy.zeros((N_total));
        particle_info['mass'] = numpy.zeros((N_total));
        particle_info['charge'] = numpy.zeros((N_total));
        particle_info['is_active_DNA'] = numpy.zeros((N_total));

        particle_info['diameter'] = numpy.zeros((N_total));
        particle_info['bonds'] = numpy.zeros((N_bonds,2));
        particle_info['bond_types'] = [];
        particle_info['bond_id'] = numpy.zeros(N_bonds);
        particle_info['running_pos_list'] = [];
        particle_info['running_diam_list'] = [];

        particle_info['Nchain_min'] = 0;
        particle_info['N_mono_min'] = 0;
        particle_info['bond_start'] = 0;

        numpy.random.seed(seed=seed_position);


        pos_params['L'] = L_A_box;
        (particle_info,position) = initialize_system_complete(N_A,seq_A,pt_info,position,particle_info,pos_params);
#        position_AW_scaffolds =

        pos_params['L'] = L;

        (particle_info,position) = initialize_system_complete(N_B,seq_B,pt_info,position,particle_info,pos_params);
        (particle_info,position) = initialize_system_complete(N_C,seq_C,pt_info,position,particle_info,pos_params);

        #Assign particle position



        #Create a snap-shot with the required particles
        snapshot = hoomd.data.make_snapshot(N=N_total, box=hoomd.data.boxdim(L=L),  particle_types=particle_types);



        snapshot.particles.position[:] = position;
        snapshot.particles.typeid[:] = particle_info['typeid'];
        snapshot.particles.mass[:] = particle_info['mass'];
        snapshot.particles.diameter[:] = particle_info['diameter'];
        snapshot.particles.charge[:] = particle_info['charge'];


        snapshot.bonds.resize(N_bonds);
        snapshot.bonds.types = particle_info['bond_types'];
        snapshot.bonds.group[:] = particle_info['bonds'];
        snapshot.bonds.typeid[:] = particle_info['bond_id'];


        # Set particle moment of inertia, diameter, and velocity

        # if solid spheres, MI approx 2/5 MR^2 for each part
        MI = numpy.multiply(2/5*snapshot.particles.mass[:],numpy.power(0.5*snapshot.particles.diameter[:],2))
        snapshot.particles.moment_inertia[0:N_total,:] =numpy.column_stack((MI,MI,MI));

        # Set thermal velocities
        V = numpy.random.normal(0.0,numpy.sqrt(T / 1.0), [snapshot.particles.N, 3]);
        Mass = numpy.column_stack((numpy.sqrt(snapshot.particles.mass[:]),numpy.sqrt(snapshot.particles.mass[:]),numpy.sqrt(snapshot.particles.mass[:])));
        snapshot.particles.velocity[:] = numpy.divide(V,Mass);

        system=hoomd.init.read_snapshot(snapshot);





        render_sphere_frame(system.take_snapshot())




        system.N_A = N_A;
        system.N_total = N_total;



        # Set up the rigid structure of the bodies
        rigid = hoomd.md.constrain.rigid();

        for ptype in particle_types:
            if is_active_DNA_list[ptype] is 1:
                
                       
                active_DNA_type = ptype;
                
                N_bs_B = int(input_params['N_bs_' + active_DNA_type +'B']);
                N_bs_C = int(input_params['N_bs_' + active_DNA_type +'C']);
                
                # Generates a function that constructs a cubic particle
                pos = generate_cubic_position(diameter_list[active_DNA_type]/3);
        
                #Generates positions of mid faces, scaled by binding site size
                face_pos = (generate_cubic_face_position(diameter_list[active_DNA_type],diameter_list[active_DNA_type+'BS']))
        
        
                pos = numpy.vstack([pos,face_pos])
        
                diam = numpy.vstack([diameter_list['AW']*numpy.ones((26,1)),diameter_list[active_DNA_type+'BS']*numpy.ones((6,1))]);
        
        
                rigid.set_param(active_DNA_type,
                                types=['AW']*26+[active_DNA_type+'BS']*N_bs_B+[active_DNA_type+'CS']*N_bs_C,
                                positions=pos[0:N_bs_B+N_bs_C+26,:],orientations=None, charges=None, diameters=diam[0:N_bs_B+N_bs_C+26]);
        

        rigid.set_param('B',
                        types=['BDNA']*1,
                        positions=[(0,0,-0.5*(diameter_list['B'])+0.5*diameter_list['BDNA'])],orientations=None, charges=None, diameters=[diameter_list['BDNA']]);

        rigid.set_param('C',
                        types=['CDNA']*1,
                        positions=[(0,0,-0.5*(diameter_list['C'])+0.5*diameter_list['CDNA'])],orientations=None, charges=None, diameters=[diameter_list['CDNA']]);



        rigid.create_bodies()








        render_sphere_frame(system.take_snapshot())




        def estimate_chain_volume(seq,diameter_list):
            volume = 0;
            for aa in seq:
                volume = volume + 4/3*numpy.pi*numpy.power(diameter_list[aa]/2,3);

            return volume;


        # Estimate the actual volume fraction

        vol_A = N_A*len(seq_A)*pow(3*diameter_list['A'],3)/float(pow(L,3));

        vol_B = N_B*(estimate_chain_volume(seq_B,diameter_list))/float(pow(L,3));
        vol_C = N_C*(estimate_chain_volume(seq_C,diameter_list))/float(pow(L,3));

        vol_fraction = vol_A+vol_B+vol_C;
        print("The volume fraction is {}".format(vol_fraction))
        print("The individual volume fraction are B = {} and C= {}".format(vol_B,vol_C))






        # Defining the inter-particle interactions
        nl = hoomd.md.nlist.cell()


        table = hoomd.md.pair.table(width=1000,nlist=nl);
        btable = hoomd.md.bond.table(width=1000);
        ci = 0;

        for typei in system.particles.types:
            cj=0;
            for typej in system.particles.types:
                table.pair_coeff.set(typei, typej, func=lj_shifted, rmin=0.1, rmax=1.25*(diameter_list[typei]+diameter_list[typej]), coeff=dict(epsilon=interaction_list[(typei,typej)], d1=diameter_list[typei],d2=diameter_list[typej],z1= charge_list[typei],z2 = charge_list[typej],lb=lb,kappa=kappa));
                cj  = cj+1;
            ci = ci+1;


        for typei in particle_info['bond_types']:
            [i1,i2] = typei.split('_');
            btable.bond_coeff.set(typei, func=harmonic, rmin=0.1, rmax=2*(diameter_list[i1]+diameter_list[i2]), coeff=dict(k=k, d1=diameter_list[i1],d2=diameter_list[i2]));



#        print(table.pair_coeff.values)




        hoomd.md.integrate.mode_standard(dt=dt,aniso=True);


        #List of all DNA center particles to not integrate
        
        DNA_core_tag_list = list(numpy.arange(N_A*len(seq_A)));
        if DNA_core_tag_list:
            print('I am here {}'.format(DNA_core_tag_list))
            groupA = hoomd.group.tag_list(name='groupA', tags=range(0,N_A*len(seq_A)));
        else:
            groupA = hoomd.group.type(name='groupA', type='A');

        #List of all rigid paricles (B and C and patches)
        rigid_all = hoomd.group.rigid();

        #List of all center particles (A,B,C) centers
        rigid_center = hoomd.group.rigid_center()

        #List of all  patches
        rigid = hoomd.group.difference(name="rigid_but_only_non-center",a=rigid_all,b=rigid_center)

        # List of all patches and A particles
        non_integrable_groups =hoomd.group.union(name="non-int",a=groupA,b=rigid);

        all = hoomd.group.all()

        # groupB = hoomd.group.type(name='groupB', type='B')
        # groupC = hoomd.group.type(name='groupC', type='C')
        #groupBC = hoomd.group.union(name="bc-particles", a=groupB, b=groupC)

        groupBC =  hoomd.group.difference(name="particles-not-typeA", a=all, b=non_integrable_groups);

        file_PE_T_KE = file_output +'_PE_T_KE.log'
        print((file_PE_T_KE))
        energy_logger = hoomd.analyze.log(filename=file_PE_T_KE,
                          quantities=['pair_table_energy', 'temperature','kinetic_energy'],
                          period=log_step,
                          overwrite=True);
        traj_dump = hoomd.dump.gsd(trajectory_file, period=log_trajectory_step, group=all, overwrite=True);



        integratorBC= hoomd.md.integrate.langevin(group=groupBC, kT=T, seed=seed);

        # This ensures that if we don't want to integrate the A particles
        # with infinite damping, then we input gamma_A <=0 for that
        if gamma_A >0:
            integratorA= hoomd.md.integrate.langevin(group=groupA, kT=T_A, seed=seed);
            integratorA.set_gamma('A',gamma_A);

        f_cluster = file_output+'_cluster_size.dat';
        logger = hoomd.analyze.log(filename=f_cluster, quantities=['cluster_size','cluster_rg','density'], period=log_step,overwrite=True)
        a = cluster_size_analyze(system)
        logger.register_callback('cluster_size', cluster_size_analyze(system))
        logger.register_callback('cluster_rg', cluster_rg_analyze(system))
        logger.register_callback('density', density_estimate(system))









        #   warm-up run if the systems are dense
        #   enough to avoid bond overlaps        
        if vol_fraction > 0.05:
           warm_up_flag =1;
           warm_up_steps =20000;
        if warm_up_flag:
            temp_energy_logger = hoomd.analyze.log(filename=file_output+'_warm_up_T.log',
                          quantities=['temperature'],
                          period=log_step,
                          overwrite=True);
            temp_energy_logger.enable();
            hoomd.md.integrate.mode_standard(dt= 0.05/k,aniso=True)
            logger.disable()
            traj_dump.disable()
            energy_logger.disable()
            T_temp= 10;
            T_threshold = 0.5;
            while (abs(T_temp-T) > T_threshold):
                hoomd.run(warm_up_steps)
                T_temp=temp_energy_logger.query('temperature');
        
            temp_energy_logger.disable();
        
        logger.enable()
        traj_dump.enable()
        energy_logger.enable()
        hoomd.md.integrate.mode_standard(dt=dt,aniso=True)



        hoomd.run(t_end);



        render_sphere_frame(system.take_snapshot())




        if drug_test:
            N_types = len(system.particles.types);
            print(N_types)
            groupCDNA = hoomd.group.type(name="cdna-particles", type='CDNA')
            for p in groupCDNA:
                if p_drug >= numpy.random.random_sample():
                    p.type = 'GS'
            groupBDNA = hoomd.group.type(name="bdna-particles", type='BDNA')
            for p in groupBDNA:
                if p_drug >= numpy.random.random_sample():
                    p.type = 'GS'



            hoomd.run(t_post_drug);


        if plot_flag:
            generate_plots(figure_output,f_cluster,file_PE_T_KE)

        if movie_flag:
            save_movie(render_sphere_frame, trajectory_file);


    print("Simulation is complete now");

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Take input and output filename and run MD simulations')
  parser.add_argument('--i',help="Name of the input_file", required = True);
  parser.add_argument('--o',help="Name of output file", required = True);
  parser.add_argument('--e',help="Name of energetics file", required = True);
  parser.add_argument('--f',help="Name of monomer identity and information file", required = True)
  parser.add_argument('--dest_folder',help="Folder to write Output data_structures too", required = False)

  parser.add_argument('--p',help="List of changing parameters", required = False)
  parser.add_argument('--outfolder',help="Name of Output folder which will be prefixed to all simulation outputs", required=False)
  args = parser.parse_args();
  # params_file = args.param_list;
  # if params_file is not None:
  #     (param_name,param_list) = params_parse(params_file)

  run_MD_simulation(args.i,args.o,args.f,args.p,args.e,args.outfolder,args.dest_folder);
