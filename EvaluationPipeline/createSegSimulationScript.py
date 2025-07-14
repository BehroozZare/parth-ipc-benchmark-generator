#Behrooz Zarebavani - 2022 - Sep
#This script has the intention to create testing suite
import os
import re
import pathlib

# IPCTestGenerator is responsible for generating batch scripts for running segmented simulation tests on a cluster.
class IPCTestGenerator:
    def __init__(self, cluster,
                 output_dir_name,
                 num_threads,
                 sing_exe,
                 sing_bind_input,
                 sing_bind_inside,
                 prog_address,
                 solver,
                 config_address,
                 simulation_dict,
                 fixed_flags,
                 test_parameters):
        """
        Initialize the test generator with cluster and simulation configuration.
        cluster: Name of the cluster (e.g., Niagara, Cedar)
        output_dir_name: Name of the output directory for scripts
        num_threads: List of thread counts to test
        sing_exe: Path to the Singularity executable
        sing_bind_input: Path to bind as input in Singularity
        sing_bind_inside: Path inside the container to bind to
        prog_address: Path to the simulation program inside the container
        solver: Solver type (e.g., MKL, CHOLMOD)
        config_address: Path to simulation config files
        simulation_dict: Dictionary of simulation parameters
        fixed_flags: Flags that are fixed for all tests
        test_parameters: Dictionary of test parameters to sweep
        """
        
        self.cluster_name = cluster
        self.output_dir = output_dir_name
        self.num_threads = num_threads
        self.sing_exe = sing_exe
        self.sing_bind_input = sing_bind_input
        self.sing_bind_inside = sing_bind_inside
        self.prog_address = prog_address
        self.solver = solver
        self.config_address = config_address
        self.sim_dict = simulation_dict
        self.fixed_flags = fixed_flags
        self.test_parameters = test_parameters
        
        self.out_path = ""
        self.time_per_sim = []
        self.verbose = True
        self.run_all_file = None
        self.test_names = []
        self.test_flags = []
        
        
    def createTests(self):
        """
        Main entry point to generate all test scripts and folders.
        """
        #Create output folder
        self.estimateTime()
        self.createOutputFolder()
        self.createRunAllScript()
        self.createTestNames()
        self.createSimulationTestScripts()
        
    
    def estimateTime(self, verbose=True):        
        """
        Estimate the required time for each simulation based on its parameters.
        """
        scaling_factor = 2.5 #For 20 cores, with trial and error it is a better number
        #Estimating time
        for sim in self.sim_dict:
            time = (self.sim_dict[sim]['Frame'] * self.sim_dict[sim]['timing'] * self.sim_dict[sim]['core'])
            h_time = int(time / 3600 / scaling_factor) + 1
            self.time_per_sim.append(h_time)
            if verbose:
                print("simulation", sim)
                print("Hours is:", self.time_per_sim[-1])
                print("*************************************************************************")


    def createOutputFolder(self):
        """
        Create the output directory structure for scripts and CSVs.
        """
        #Building the directories and add scripts
        parent_dir = os.getcwd()
        self.out_path = os.path.join(parent_dir, 'scripts', self.output_dir)
        pathlib.Path(self.out_path).mkdir(parents=True, exist_ok=True) 
        csv_dir = os.path.join(self.out_path, "csv")
        pathlib.Path(csv_dir).mkdir(parents=True, exist_ok=True) 
        
        
    def createRunAllScript(self):
        """
        Create a master script to submit all generated batch jobs.
        """
        #Create run all script
        self.run_all_file = open(os.path.join(self.out_path, 'run_all.sh'), "w")
        self.run_all_file.write("#!/bin/bash\n\n\n\n")        
        
    def createTestNames(self):
        """
        Generate all combinations of test names and flags based on parameters.
        """
        #Create initial test name and flags
        test_name_prev = []
        temp_flag_prev = []
        for th in self.num_threads:
            test_name_prev.append("numThreads_" + th + "_" + "SolverType_" + self.solver)
            temp_flag_prev.append(self.fixed_flags  + " --SolverType=" + self.solver + " --numThreads=" + th)
        
        
        self.test_names = []
        self.test_flags = []
        
        #Creating the names
        for par_name in self.test_parameters:
            par_values = self.test_parameters[par_name]
            self.test_names = []
            for test in test_name_prev:
                for value in par_values:
                    self.test_names.append(test + "_" + par_name + "_" + value)
            test_name_prev = self.test_names
            
        #Creating the flags
        for par_name in self.test_parameters:
            par_values = self.test_parameters[par_name]
            self.test_flags = []
            for test in temp_flag_prev:
                for value in par_values:
                    self.test_flags.append(test + " --" + par_name + "=" + value)
            temp_flag_prev = self.test_flags
        
        print(len(self.test_names), len(self.test_names))
    def getParametersFromName(self, name):
        """
        Parse a test name string into a dictionary of parameter values.
        """
        flag_dict = {}
        temp = name.split("_")
        for cnt in range(int(len(temp) / 2)):
            flag_dict[temp[2 * cnt]] = temp[2 * cnt + 1]
        return flag_dict
                
    
    def createSimulationTestScripts(self):
        """
        For each simulation and test, generate a batch script to run the simulation with the correct parameters.
        """
        sim_cnt = -1
        for sim in self.sim_dict:
            sim_cnt = sim_cnt + 1
            sim_path = os.path.join(self.out_path, sim)
            pathlib.Path(sim_path).mkdir(parents=True, exist_ok=True) 
        
            for cnt in range(len(self.test_names)):  
                test = self.test_names[cnt]
                flag_dict = self.getParametersFromName(test)
                run_method_script = open(os.path.join(sim_path, test + '.sh'), "w")
                #adding 
                run_method_script.write("#!/bin/bash\n\n\n\n")

                run_method_script.write("#SBATCH --account=rrg-mmehride\n")
                run_method_script.write('#SBATCH --job-name="IPC_'+ sim + "_" + test + '"\n')
                run_method_script.write('#SBATCH --output="IPC_'+ sim + "_" + test + '_%j.out"\n')
                run_method_script.write('#SBATCH --error="IPC_'+ sim + "_" + test + '_%j.error"\n')
                run_method_script.write("#SBATCH --nodes=1\n")
                if self.cluster_name=="Niagara":
                    run_method_script.write("#SBATCH --ntasks-per-node=40\n")
                    run_method_script.write("#SBATCH --export=ALL\n")
                    run_method_script.write("#SBATCH -t " + str(min(self.time_per_sim[sim_cnt],24)) + ":00:00\n")
                    run_method_script.write("#SBATCH --constraint=cascade\n\n\n\n")
                elif self.cluster_name == "Cedar":
                    run_method_script.write("#SBATCH --ntasks-per-node=48\n")
                    run_method_script.write("#SBATCH --export=ALL\n")
                    run_method_script.write("#SBATCH -t " + str(min(self.time_per_sim[sim_cnt],24)) + ":00:00\n")
                    run_method_script.write("#SBATCH --constraint=cascade\n\n\n\n")
                else:
                    print("Please implement the batch config for cluster", self.cluster_name)
                
                run_method_script.write('export SINGULARITY_BIND="' + self.sing_bind_input + ':/' + self.sing_bind_inside + '"\n')
                run_method_script.write('export SING_SIF=/"' + self.sing_exe + '"\n')
                run_method_script.write('export num_threads=' + flag_dict["numThreads"] + "\n")
                run_method_script.write('export MKL_NUM_THREADS=$num_threads\n')
                run_method_script.write('export OMP_NUM_THREADS=$num_threads\n')
                run_method_script.write('export VECLIB_MAXIMUM_THREADS=$num_threads\n')
                run_method_script.write('export PROG_PATH=/' + self.sing_bind_inside + self.prog_address + '\n')
                run_method_script.write('export Config=/' + self.sing_bind_inside + self.config_address + sim + ".txt\n")
                run_method_script.write('export progMode=100\n\n\n')

                run_method_script.write('singularity exec $SING_SIF bash -c "$PROG_PATH $progMode $Config' + " --SimName=../../csv/" + sim + "_" + test + ' --output=/mnt/' + os.path.join(self.output_dir, sim, test) + self.test_flags[cnt] + '"\n')

                run_method_script.close()
                self.run_all_file.write('sbatch  ' + os.path.join(sim, test + '.sh') + '\n')
            self.run_all_file.write('\n\n\n')
        self.run_all_file.close()


def createSegConfig(status_file_path, input_address, output_address, sim_name, total_frame, segments):
    """
    Generate segmented config files for a simulation, splitting the total frames into segments.
    Each segment gets its own config file with updated time and restart information.
    """
    #read config file change the time and add the restart and generate the config
    cnt = 0
    for seg in range(segments):    
        output_config_file = open(os.path.join(output_address, sim + '_seg' + str(cnt) + '.txt'), "w")
        input_config_file = open(os.path.join(input_address, sim + '.txt'), "r")
    
        start_frame = int(cnt * (total_frame / segments))
        end_frame = int((cnt + 1) * (total_frame / segments))
        cnt = cnt + 1
        num_frame = int(end_frame) - int(start_frame)
        for line in input_config_file:
            if line.startswith('time '):
                stepsize = float(line.split(" ")[-1])
                output_config_file.write('time ' + str((int(start_frame) + num_frame + 1) * stepsize) + ' ' + str(stepsize) + '\n')
            else:
                output_config_file.write(line)
        output_config_file.write('\nrestart' + ' ' + status_file_path + sim + '/status' + str(start_frame) + "\n")
        
        output_config_file.close()
        input_config_file.close()
    pass

# The simulations dictionary defines the set of simulations and their resource/time requirements.
simulations = {
    '8_rollerBall': {'Frame': 1000, 'memory': 357, 'timing':  200, 'Iteration': 39.7, 'core': 4},
    '13_dolphinFunnel': {'Frame': 800, 'memory': 357, 'timing':  100, 'Iteration': 39.7, 'core': 4},
    '4_rodsTwist': {'Frame': 4000, 'memory': 2638, 'timing':  141.5, 'Iteration': 2.8, 'core': 8},
    '1_squeezeOut': {'Frame': 1500, 'memory': 1700, 'timing':  252, 'Iteration': 42.5, 'core': 8},
    '16_armaRoller_E1e5': {'Frame': 400, 'memory': 3651, 'timing':  346, 'Iteration': 66.8, 'core': 4},
    '12_matOnBoard': {'Frame': 200, 'memory': 357, 'timing':  100, 'Iteration': 39.7, 'core': 4},
    '14_matTwist': {'Frame': 2500, 'memory': 4546, 'timing':  776.2, 'Iteration': 34.5, 'core': 8},
              }

if __name__ == "__main__":
    # Main script entry: generate segmented config files and test scripts for all simulations and solvers.
    NUM_SEG = 10
    for sim in simulations:
        createSegConfig(
            status_file_path = "/mnt/Checkpoints/",
            input_address = "/home/behrooz/Desktop/IPC_Project/IPC-dev2/input/paperExamples/",
            output_address = "/home/behrooz/Desktop/IPC_Project/IPC-dev2/input/seg/",
            sim_name = sim,
            total_frame = simulations[sim]['Frame'],
            segments=NUM_SEG)
                        
    seg_simulations = {}
    for sim in simulations:
        for i in range(NUM_SEG):
            seg_simulations[sim + '_seg' + str(i)] = {'Frame': simulations[sim]['Frame'] / NUM_SEG, 'memory': 1700, 'timing':  simulations[sim]['timing'], 'Iteration': simulations[sim]['Iteration'], 'core': 8}
    
    
    strum_tester = IPCTestGenerator(cluster="Niagara",
                output_dir_name="MKL_SIM_Test",
                num_threads=["20"],
                sing_exe="$SCRATCH/DOCKER/PARTH_DOCKER.sif",
                sing_bind_input="$SCRATCH/DOCKER/Source",
                sing_bind_inside="mnt",
                prog_address = "/IPC-dev2/build/IPC_bin",
                solver="MKL",
                config_address="/IPC-dev2/input/seg/",
                simulation_dict=seg_simulations,
                fixed_flags=" --DoAnalysis=0",
                test_parameters={"IM":['0']})
    strum_tester.createTests()
    
    
    strum_tester = IPCTestGenerator(cluster="Niagara",
                output_dir_name="EIGEN_SIM_Test",
                num_threads=["20"],
                sing_exe="$SCRATCH/DOCKER/PARTH_DOCKER.sif",
                sing_bind_input="$SCRATCH/DOCKER/Source",
                sing_bind_inside="mnt",
                prog_address = "/IPC-dev2/build/IPC_bin",
                solver="EIGEN",
                config_address="/IPC-dev2/input/seg/",
                simulation_dict=seg_simulations,
                fixed_flags=" --DoAnalysis=0",
                test_parameters={"IM":['0']})
    strum_tester.createTests()
    
    
    
    strum_tester = IPCTestGenerator(cluster="Niagara",
                output_dir_name="PARSY_SIM_Test",
                num_threads=["20"],
                sing_exe="$SCRATCH/DOCKER/PARTH_DOCKER.sif",
                sing_bind_input="$SCRATCH/DOCKER/Source",
                sing_bind_inside="mnt",
                prog_address = "/IPC-dev2/build/IPC_bin",
                solver="PARSY",
                config_address="/IPC-dev2/input/seg/",
                simulation_dict=seg_simulations,
                fixed_flags=" --DoAnalysis=0",
                test_parameters={"IM":['0']})
    strum_tester.createTests()
    
    

    strum_tester = IPCTestGenerator(cluster="Niagara",
                    output_dir_name="CHOLMOD_CHECKPOINT_Test",
                    num_threads=["20"],
                    sing_exe="$SCRATCH/DOCKER/PARTH_DOCKER.sif",
                    sing_bind_input="$SCRATCH/DOCKER/Source",
                    sing_bind_inside="mnt",
                    prog_address = "/IPC-dev2/build/IPC_bin",
                    solver="CHOLMOD",
                    config_address="/IPC-dev2/input/seg/",
                    simulation_dict=seg_simulations,
                    fixed_flags=" --DoAnalysis=0",
                    test_parameters={"IM":['0']})
    strum_tester.createTests()
    