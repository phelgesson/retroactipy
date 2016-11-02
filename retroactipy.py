# Built-in modules #
import os, re, sys
from datetime import datetime
#
# Third party modules #
import numpy as np
import sh, shutil
from subprocess import call
#
# Personal modules
sys.path.append("/home/petter/pythonWork") 
import plotter as pr # only use for pr.show() within project
#


# Short forms #
npcat = np.concatenate
#

class Nuclide(object):
    """Class describing the basic features of a nuclide"""
    
    def __init__(self, atomic_number, element, mass, state):
        """Constructor of Nuclide objects"""
        self.atomic_number = int(atomic_number)
        self.element = element
        self.mass = int(mass)
        self.state = int(state)
        
        

class CrossSection(object):
    """"""
    
    def __init__(self, mt, energies, cross_section):
        """Constructor of CrossSection objects"""
        self.mt = int(mt)
        mts = np.array([1, 2, 102, 18])
        xs_types = ['total', 'elastic', 'capture', 'fission']
        xs_strings = ['(n,tot)', '(n,el)', '(n,$\gamma$)', '(n,f)']
        idx = np.where(mt == mts)[0]
        self.sammy_type = xs_types[idx]
        self.string = xs_strings[idx]
        self.energies = np.array(energies)
        self.cross_section = np.array(cross_section)
        self.experiment = None
        
        
class PseudoExperiment(object):
    
    
    def __init__(self, cross_sections):
        self.cross_sections = cross_sections
        for xs in cross_sections:
            xs.experiment = self
        self.generate_cov()
        
    def __getitem__(self, idx):
        """Uses mt number or type for indexing."""
        if type(idx) == list:
            answer = []
            for i in idx:
                answer.append(self[i])
            return answer
        else:
            for xs in self.cross_sections:
                if xs.mt == idx :
                    return xs
            raise IndexError(
                'mt or type = %s not found in PsedoExperiment instance.' 
                % str(idx))
        
    def generate_cov(self):
        for xs in self.cross_sections:
            xs.cov = np.diag(xs.cross_section*.10) # dummy
        
        
        
class EndfFile(object):
    
    def __init__(self, path):
        """Constructor with error checks for the path"""
        path = absolute_path(path)
        assert path != ''
        if not os.path.exists(path):
            raise ValueError('Path %s does not exist' % path)
        assert not os.path.isdir(path)
        self.path = path

        
class CompleteEndfFile(EndfFile):
    """"""
    
    def __init__(self, path):
        EndfFile.__init__(self,path)
    
    def create_resonances(self, path, pathcov): 
        """creates ResonanceFile instance and produces the 
        corresponding file(s) containing the subsection of the 
        ENDF file with resonance parameters and (if available) 
        covariances"""
        try: 
            self.nuclide 
        except AttributeError:
            self.nuclide = self.create_nuclide()
        f = open(self.path,'r')
        fres = open(path,'w')
        frescov = open(pathcov,'w')
        try:
            found_res = False
            found_cov = False
            while True:
                line = f.readline()
                if line == '': break
                el = EndfLineParser(line)
                # Look for MF = 2 and MF = 32
                if el['mf'] == 2 and el['mt'] == 151:
                    # Write out
                    if not found_res:
                        fres.write(old_line) # starts one line above
                        found_res = True
                    while el['mf'] == 2 and el['mt'] == 151:
                        fres.write(el.line)
                        el = EndfLineParser(f.readline())
                    fres.write(el.line) # two more
                    fres.write(f.readline()) 
                    #
                    # Create corresp. ResonanceFile object
                    res_file = ResonanceFile(path, self.nuclide)
                elif el['mf'] == 32 and el['mt'] == 151:
                    # Write out
                    if not found_cov:
                        fres.write(old_line)
                        found_cov = True
                    while el['mf'] == 32 and el['mt'] == 151:
                        fres.write(el.line)
                        el = EndfLineParser(f.readline())
                    fres.write(el.line) # two more
                    fres.write(f.readline()) 
                    # Create corresp. ResonanceCovFile object
                    res_file.cov = ResonanceCovFile(pathcov)
                elif found_res and found_cov:
                    break # no need to continue
                #
                old_line = line
            #If no MF = 2 found, raise exception
            if not found_res:
                os.remove(path)
                raise ValueError('No resonance "file" in %s' % self.path)
            if not found_cov:
                os.remove(pathcov)
        finally:
            f.close()
            fres.close()
            frescov.close()
        return res_file
            
                        
    def create_nuclide(self):
        """Parses head of the file and calls create resonances,
        returns Nuclide instance"""
        with open(self.path) as f:
            f.readline() # Skip zeroth line
            # line 1
            el = EndfLineParser(f.readline())             
            za = int(el[0]) # float combining atomic number and mass
            neutron_mass = 1.00866491588
            weight = el[1] * neutron_mass
            atomic_number = za / 1000
            mass = za - atomic_number * 1000
            mat_nbr = el['mat']
            #
            # line 2
            el = EndfLineParser(f.readline())
            state = el[2]
            #
            for k in range(3,6): # go to 5th line
                el = EndfLineParser(f.readline())
            element = el[0].split('-')[1]
            self.nuclide = Nuclide(atomic_number, element, mass, state)
            self.nuclide.mat_nbr = mat_nbr
            self.nuclide.weight = weight
            self.nuclide.channel_radius = 0.
            return self.nuclide
     
        
class ResonanceFile(EndfFile):
    
    def __init__(self, path, nuclide):
        EndfFile.__init__(self, path)
        self.nuclide = nuclide
        
        
        
class ResonanceCovFile(EndfFile):
    
    def __init__(self, path):
        EndfFile.__init__(self, path)
        
        
class Main(object):
    """God."""
    
    def __init__(self, endf_path, directory = None, sammy_executable = 'sammy'):
        """"""
        # Error checks #
        endf_path = absolute_path(endf_path)
        assert os.path.exists(endf_path)
        assert not os.path.isdir(endf_path)
        if sh.which(sammy_executable) == None:
            raise ValueError("Command '%s' not found." % sammy_executable +
            "Please set attribute 'sammy_executable'.")
        #
        self.endf_file = CompleteEndfFile(endf_path)
        # If directory == None, set to the directory of the endf file #
        directory = absolute_path(directory)
        if directory == None:           
            directory = re.sub('/[^/]+$', '', path)
        else:
            os.makedirs(directory)
        #
        respath = directory + 'mf2' # change
        rescovpath = directory + 'mf32' # change
        self.nuclide = self.endf_file.create_nuclide()
        
        
    def add_experiment(self, energies):
        self.experiment = self.nuclide[0].reconstruct_cross_sections(energies)
        



class EndfLineParser(object):
    """
    Class which is initialized using a line from an ENDF file, and then can be
    used to call for certain parts of that line, using square brackets. 
    Indices of normally formatted lines (or oddly formatted using self.splits)
    can be used, or certain keywords. Lists combining these two, and slices, are 
    also allowed.
    """
    
    def __init__(self, line, splits = None):
        self.line = line
        self.default_splits = npcat((np.array(range(0,7,1))*11,[70,72,75,80]))
        if splits == None:
            self.splits = self.default_splits
        
    def __getitem__(self,idx):
        tp = type(idx)
        if tp == int:
            idx = idx % len(self)
            word = self.line[self.splits[idx]:self.splits[idx+1]]
            if word.strip().isdigit():
                # Integer if "possible"
                return int(word)
            try:
                # Try to return float #
                # identify +/- signs that are powers: #
                m = re.findall('([0-9.]+)([+-][0-9]+)',word)
                if m != []: # power found
                    # reconstruct with E for power:
                    return float(m[0][0] + 'E' + m[0][1])
                #
                elif re.search('[^ \t\n\r\f\v]',word) == None:
                    return 0. # empty entry => 0
                else:
                    return float(word)
                #
            except ValueError:
                # if not a number, return string
                return word
        elif tp == list:
            answer = []
            for i in idx:
                answer.append(self[i])
            return answer
        elif tp == slice:
            if idx.start == None:
                start = 0
            else:
                start = idx.start
            if idx.stop == None:
                stop = len(self)
            else:
                stop = idx.stop
            if idx.step == None:
                step = 1
            else:
                step = idx.step
            return self[range(start, stop, step)]
        elif tp == str:
            self.splits = self.default_splits
            idx = idx.lower()
            if idx in ['mf', 'file']:
                return self[-3]
            elif idx in ['mt', 'tape']:
                return self[-2]
            elif idx[:3] == 'mat':
                return self[-4]
            elif idx[:3] in ['num', 'nbr']:
                return self[-1]
            else:
                raise KeyError('Invalid key: %s' % idx)
        else:
            raise TypeError('%s invalid as index.' % str(tp))
        
        
    def __len__(self):
        return len(self.splits)-1
        
'''
class EndfParser(object):
    
    def __init__(self,file_handle):
        self.file_handle = file_handle
        
    def __getitem__(self):
        yield EndfLineParser(self.file_handle.readline())
        
'''

class SammyRunner(object):
    """Class collecting methods and settings (attributes) for running SAMMY"""
    
    def __init__(self, 
                executable = 'sammy', 
                data_format = 'twenty',
                temperature = 0., # This should be changed in future work
                default_thickness = 0.003,
                baby_steps = True,
                cleanup = True,
                convergence_dict = dict(
                                        max_iter = 200,
                                        length = 10, 
                                        tol = 0.01, # reduced chi-sq tol
                ),
                lcomp = 2,
                ndigits_f32 = 6,
        ):
        
        # Settings are stored as attributes #
        self.executable = executable
        self.data_format = data_format
        self.temperature = temperature
        self.default_thickness = default_thickness
        self.default_xs_type = 'total'
        self.cleanup = cleanup
        self.convergence_dict = convergence_dict
        self.convergence_dict['min_iter'] = self.convergence_dict['length']
        self.baby_steps = baby_steps
        self.lcomp = lcomp
        self.ndigits_f32 = ndigits_f32
        #
        
        
    def endf2inp_par_ndf(self, resonance_file, out_paths, xs_measurement = None,
                     flag_all = False):
        """
        Given a ResonanceFile instance and a list with the paths to desired .inp
        .res and .ndf files, the function uses SAMMY to produce these files.
        """
        nucl = resonance_file.nuclide
        for k in range(len(out_paths)):
            out_paths[k] = absolute_path(out_paths[k])
        # name of temporary files
        tempinp = temp_file_gen('Sammy_endf2par', 'inp')
        tempndf = temp_file_gen('Sammy_endf2par', 'ndf') 
        # 
        keylines = [
                    'do not solve bayes equations',
                    'input is ENDF/B file MAT=%i' % nucl.mat_nbr,
                    'use energy range from ENDF/B file 2',
        ]
        if len(out_paths) > 2: 
            keylines.insert(0,'automatic ndf file creation')
            keylines.append('ndf file is in key word format')
        title = 'SAMMY inp (generated from) for ENDF '
        if len(title) + len(resonance_file.path) > 80: 
            pathstr = '...' + resonance_file.path[-(77 - len(title)):]
        else: pathstr = resonance_file.path
        title += pathstr
        self.construct_inp(tempinp, nucl, keylines, title, xs_measurement)
        desired = ['SAMNDF.INP', 'SAMNDF.PAR', 'SAMNDF.NDF'][:len(out_paths)]
        self.run([tempinp, resonance_file.path, 'dummy.dat']
        , desired, out_paths)
        
        
        # Modify input for use in other steps
        keyremove = ['do not solve bayes equations',
                     'automatic ndf file creation']
        keylines = []
        if flag_all: keylines.append('flag all resonance parameters')
        self.modify_inp(out_paths[0], keylines, keyremove = keyremove)
        #
        # Change LCOMP to self.lcomp and set NDIGITS to self.ndigits_f32
        if len(out_paths) > 2:
            with open(out_paths[2],'r') as f1:
                with open(tempndf,'w') as f2:
                    for line in f1:
                        if line[:5].lower() == 'lcomp':
                            f2.write('Lcomp = %i\n' % self.lcomp)
                            if self.lcomp == 2:
                                f2.write('Ndigits = %i\n' % self.ndigits_f32)
                        else:
                            f2.write(line)
        shutil.move(tempndf, out_paths[2])
        #
        # Clean up #
        if self.cleanup:
            os.remove(tempinp)
        #
        
        
    def inp_par_ndf_cov2endfs(self, input_paths, outpaths):
        """
        File 32 'alone is not capable...' See SAMMY-8 manual rev 7 p. 532
        """
        keylines = [
                    'do not solve bayes equations',
                    'ENDF/B-VI file 2 is wanted',
                    'ndf file is in key word format',
                    'put covariance matrix into endf file 32'
                   ]
        
        
        tempinp = temp_file_gen('Sammy_inps2endfs')
        shutil.copy(input_paths[0], tempinp)
        input_paths[0] = tempinp
        self.modify_inp(input_paths[0], keylines)
        self.run(input_paths, ['SAMMY.NDF', 'SAMMY.N32'], outpaths)
        if self.cleanup:
            os.remove(tempinp)
        
        
    def reconstruct(self, 
    				resonance_file, 
		    		energies # np.array with energies, or list of such arrays
		    		):
        """
        Reconstructs the cross sections using parameters in resonance_file 
        object at the provided energies, and using the 'settings' of the 
        SammyRunner object. As of now, the reconstructions are done at SAMMYs 
        default energies and interpolated, this should be changed and should 
        enable Doppler Broadening (self.temperature > 0.).
        """
        # Local variables, not prone to change #
        mt_order_in_lst = [1,2,102]
        #
        # Error checks #
        if self.temperature > 0.: 
            raise ValueError('Doppler broadening not available yet.')
        try: energies = np.array(energies) 
        except: raise ValueError('Invalid energies.')
        #
        # Set temporary file names #
        tempbase = temp_file_gen('Sammy_reconstruct')
        tempinp = tempbase + '.inp'
        temppar = tempbase + '.par'
        templst = tempbase + '.lst'
        #
        # Generate input files from ENDF
        self.endf2inp_par_ndf(resonance_file, [tempinp, temppar])
        #
        # Prepare actual reconstruction run by modifying input #
        newkeylines = ['do not solve bayes equations'
                       'reconstruct cross sections',
                       'ignore input binary covariance file',
            ]
        self.modify_inp(tempinp, newkeylines)
        #
        # Run SAMMY to reconstruct  cross sections #
        self.run([tempinp,temppar,'dummy.dat'],['SAMMY.LST'],[templst])
        #
        # Read output, generate list of CrossSection instances #
        cross_section_array = self.read_data(templst)
        cross_sections = []
        E = energies
        for k in range(len(mt_order_in_lst)):
            if type(energies) == list and len(energies) == len(mt_order_in_lst):
        		E = energies[k]
            cross_sections.append(CrossSection(mt_order_in_lst[k], E, np.interp(
            	E, cross_section_array[:,0], cross_section_array[:,k+1])))
        #
        # Clean up #
        if self.cleanup: 
            for p in [tempinp, temppar, templst]: os.remove(p)
        #    
        return cross_sections
        
    def fit(self, resonance_file, experiment, outpaths):
        """
        Produces SAMMY input and runs SAMMY 'without prior',
        returns a new ResonanceFile instance with resonance parameters fit 
        to the pseudo-experiment. The fit is so far simplified such that 
        'true GLSQ' is used for (n,tot) and then Bayesian updating is used 
        for the other channels. This would be equivalent to GLSQ all the way   
        if there were no non-linearities and no correlations between the
        different reactions. Another very important simplification: the exp.
        data can only include uncorrelated uncertainties, so far.
        All this has to be developed further later (outside course).
        """
        
        inp = temp_file_gen('Sammy_fit','inp')
        par = temp_file_gen('Sammy_fit','par')
        cov = temp_file_gen('Sammy_fit','cov')
        ndf = temp_file_gen('Sammy_fit','ndf')
        
        parout = temp_file_gen('Sammy_fit','out.par')
        covout = temp_file_gen('Sammy_fit','out.cov')
        
        self.endf2inp_par_ndf(resonance_file, [inp, par, ndf], 
            experiment[1], flag_all = True)
        self.modify_inp(inp, keyremove = ['mlbw formalism is wanted'])
        
        self.g_least_squares(inp, par, experiment[1],
            parout, covout)
        shutil.move(parout, par)
        shutil.move(covout, cov)
        
        
        for xs in experiment[[2, 102]]:
            self.bayesian([inp, par, cov], xs, 
                          [parout, covout])
            shutil.copy(parout, par)
            shutil.copy(covout, cov)
        
        self.inp_par_ndf_cov2endfs([inp, par, ndf, cov], outpaths)
        
        if self.cleanup:
            for p in [inp, par, cov, ndf]: os.remove(p)
        
        resonance_file_out = ResonanceFile(outpaths[0], resonance_file.nuclide)
        resonance_file_out.cov = ResonanceCovFile(outpaths[1])
        return resonance_file_out
        
    def g_least_squares(self, 
             sammy_inp, sammy_par, # assumed out_paths of self.endf2inp_par_ndf()
             xs_measurement, # PseudoExperiment instance
             out_par, out_cov
             ):
        """
        Performs GLSQ fit to one type of cross section using SAMMY itertively 
        until convergence. 
        """
        # Set up out paths
        final_outs = [out_par, out_cov]
        for k in range(len(final_outs)):
            final_outs[k] = absolute_path(final_outs[k])
        #
        
        # Set up and move to temporary directory #
        inp_par = [sammy_inp, sammy_par]
        originals = ['original.' + s for s in ['inp', 'par']]
        directory = temp_file_gen('Sammy_glsq', directory = True)
        origdir = os.getcwd()
        for k in range(len(originals)):    
            shutil.copy(inp_par[k],directory + originals[k])
        os.chdir(directory)
        #
        try: # Combined with finally to ensure moving back to original directory
            # Set file names #
            dat_path =  'data.dat'
            inp_path = 'inp.inp'
            #
            # Construct .dat file #
            data = np.array([xs_measurement.energies, 
                            xs_measurement.cross_section,
                            np.diag(xs_measurement.cov)]).T
            self.write_data(data, dat_path)
            #
            # Modify input 
            shutil.copy(originals[0], inp_path)
            keylines = [
                        'use least squares to give pcm for initial parameters',
                        'chi squared is wanted',
                        'ignore input binary covariance file',
                        ]
            keyremove = ['do not solve bayes']
            if self.baby_steps:
                keylines.append('take baby steps with least squares method')
            self.modify_inp(inp_path, keylines, keyremove,
                new_type = xs_measurement.sammy_type)
            #
            # Iteratively run SAMMY with this .inp, .dat and updating .par,
            # until convergence or self.convergence_dict['max_iter'] is hit
            inp_paths = [inp_path, '', dat_path]
            desired_outputs = ['SAMMY.PAR', 'SAMMY.LPT', 'SAMMY.COV']
            chisq = []
            convdict = self.convergence_dict
            for k in range(convdict['max_iter']):
                # Set file paths #
                current_nbr_str = str(k).zfill(3)
                current_par = 'par%s.par' % current_nbr_str
                output = 'out%s.lpt' % current_nbr_str
                cov = 'out.cov'
                next_par = 'par%s.par' % str(k + 1).zfill(3)
                #
                # Set input to self.run()
                output_names = [next_par, output, cov]
                inp_paths[1] = current_par
                if k == 0:
                    shutil.copy(originals[1], current_par)
                #
                # Run SAMMY
                self.run(inp_paths, desired_outputs, output_names)
                #
                # Find number of points and chi-squared in output
                with open(output,'r') as f: 
                    for line in f:
                        if "Number of experimental data points =" in line:
                            npoints = int(line.split("=")[1])
                        elif "CUSTOMARY CHI SQUARED =" in line:
                            chisq.append(float(line.split("=")[1]))
                            break
                    else:
                        raise ValueError('No chi^2 found in SAMMY output.')
                #
                # Check if converged
                if k > max(convdict['min_iter'], convdict['length']):
                    last_chisqs = np.array(chisq[-convdict['length']:]
                        ).reshape((convdict['length'],1))
                    max_diff = np.max(np.abs(last_chisqs - last_chisqs.T))
                    if max_diff < convdict['tol']: # convergence reached
                        message = 'Converged after %i iterations' % (k + 1)
                        break
                #
            else: # reached if no break statement
                message = 'Did not converge.'
            #
            # Treat converged case
            local_outs = [next_par, cov]
            self.write_data(np.array(chisq), 'chisq.txt')
            for j in range(len(local_outs)):
                shutil.copy(local_outs[j], final_outs[j])
            #
        finally:
            os.chdir(origdir)
            
        # Clean up #
        if self.cleanup: sh.rm('-rf', directory)
        #
        return message
        
    def bayesian(self, 
                 sammy_inputs, # inp (from self.endf2inp_par_ndf()), par, cov
                 xs_measurement, # CrossSection instance
                 sammy_outputs # par, cov
                ):
        """
        Performs one run of SAMMY, performing Bayesian update of resonance
        parameters in sammy_par.
        """
        tempinp = temp_file_gen('Sammy_bayesian','inp')
        tempdata = temp_file_gen('Sammy_bayesian','dat')
        
        shutil.copy(sammy_inputs[0], tempinp)
        sammy_inputs[0] = tempinp
        sammy_inputs.insert(2,tempdata)
        
        data = np.array([xs_measurement.energies, 
                        xs_measurement.cross_section,
                        np.diag(xs_measurement.cov)]).T
        self.write_data(data, tempdata)
        
        
        self.modify_inp(sammy_inputs[0], new_type = xs_measurement.sammy_type)
        self.run(sammy_inputs, ['SAMMY.PAR','SAMMY.COV'], 
                 output_names = sammy_outputs)
        
        if self.cleanup:
            for p in [tempinp, tempdata]: os.remove(p)
        
        
            
    def construct_inp(self, path, nuclide, keylines, title = None, 
            xs_measurement = None):
        """
        Given certain information, in 'self', 'nuclide', and 'keylines',
        the method constructs an input file at 'path'.
        """
        keylines = [self.data_format] + keylines 
        try: thickness = xs_measurement.experiment.thickness 
        except: thickness = self.default_thickness
        try: xs_type = xs_measurement.sammy_type
        except: xs_type = self.default_xs_type
        if title == None:
            title = 'SAMMY inp file for %s-%s' % (nuclide.element,nuclide.mass)
        with open(path,'w') as f:
            f.write(title[:80] + '\n')
            write_to_position(f, ['%s%i' % (nuclide.element, nuclide.mass), 
                nuclide.weight], [0, 10])
            for keyline in keylines + ['']:
                f.write(keyline + '\n')
            write_to_position(f, self.temperature, 0)
            write_to_position(f, [nuclide.channel_radius, thickness],[0, 10])
            write_to_position(f, xs_type, 0)

    def modify_inp(self, path, keylines = [], keyremove = [], new_type = None):
        tempnew = temp_file_gen('Sammy_modify_inp')
        keylines_complete = False
        with open(path,'r') as f1:
            with open(tempnew,'w') as f2:
                for line in f1:
                    if line.strip() == '':
                        for keyline in keylines:
                            f2.write(keyline + '\n')
                        keylines_complete = True
                    if (new_type != None and line[:3].lower() in 
                               ['cap','tot','ela','tra'] and keylines_complete):
                        f2.write(new_type + '\n')
                    else:
                        if not any([line.lower()[:len(keyrem)] == keyrem.lower()
                                for keyrem in keyremove + keylines]):
                            f2.write(line)
        shutil.move(tempnew, path)
            
    def run(self, input_paths, desired_outputs, output_names = None):
        """
        Given a list of paths to input files (posibly dummies), desired output
        files and desired names of these output files, the method constructs
        a shell script which runs SAMMY and answers the prompt such that it is
        fed with the input, and then obtains the desired output files at their 
        respective paths.
        """
        tempdir = temp_file_gen('Sammy_run', directory = True)
        if output_names == None: output_names = desired_outputs
        assert len(desired_outputs) == len(output_names)
        for k in range(len(input_paths)):
            # SAMMY doesn't accept too long variable names #
            try: shutil.copy(input_paths[k],tempdir)
            except IOError: pass # may be dummy file
            input_paths[k] = re.sub('.*/', '', input_paths[k])
            # 
        for k in range(len(output_names)):
            output_names[k] = absolute_path(output_names[k])
        origdir = os.getcwd()
        os.chdir(tempdir)
        try: # to ensure moving back to original directory
            run_file = 'run_file'
            with open(run_file, 'w') as f:
                f.write('#! /bin/sh\n\nsammy<<EOF\n')
                for path in input_paths + ['','EOF']:
                    f.write(path + '\n')
            sh.chmod('755', run_file)
            call('./%s > terminal.out 2>&1' % run_file, shell = True)
            for k in range(len(desired_outputs)):
                shutil.move(desired_outputs[k],output_names[k])
        finally:
            os.chdir(origdir)
            
        if self.cleanup:
            sh.rm('-rf',tempdir)
            
    def read_data(self, path):
        """Reads data in 'path' using format in 'self.data_format', returns 
           numpy array."""
        if self.data_format == 'twenty': 
            length = 20
        else: raise ValueError("self.data_format = '%s' unknown." % 
            self.data_format)
        data = []
        with open(path,'r') as f:
            for line in f:
                data.append([float(line[k:(k + length)]) for k in range(
                    0, len(line.strip('\n')),length)])
        return np.array(data)
    
    def write_data(self, data, path):
        """Writes numpy array 'data' into 'path' using format in 
           'self.data_format'."""
        if self.data_format == 'twenty': 
            length = 20
        else: raise ValueError("self.data_format = '%s' unknown." % 
            self.data_format)
        if len(data.shape) == 1: data = data.reshape((data.shape[0],1))
        with open(path,'w') as f:
            for k in range(data.shape[0]):
                f.write(''.join(
                    [str(data[k,l]).rjust(length) for l in range(data.shape[1])]
                    ) + '\n')


### Misc. functions ###        
def absolute_path(path):
    """
    Returns the absolute path given a relative (or absolute) path. Also handles
    the use of '~' for the home directory.
    """
    path = re.sub('~', os.environ['HOME'], str(path))
    if path[0] != '/':
        path = str(sh.pwd()).strip() + '/' + path
    return path
    
def temp_file_gen(base, extension = '', directory = False, always_first = 'temp_', time = False):
    """Generates name of temporary file including a number,
    creates empty file."""
    base = re.sub('[^/]*$','',base) + always_first + re.sub('.*/','',base).strip('_')
    extension = '.'[:(len(extension)*(extension[:1] != '.'))] + extension
    k = 0
    while True:
        if time: tstr = '_' + datetime.now().strftime('%m%d%H%M%S')
        else: tstr = ''
            
        attempt = base + tstr + '_' + str(k).zfill(3) + extension
        if not os.path.exists(attempt):
            if directory:
                os.makedirs(attempt)
                attempt += '/'
            else:
                open(attempt,'w').close()
            return attempt
        k += 1
        
def write_to_position(f, # handle of open file
                      words, # any object with a string method, or list of such
                      positions, # integer or list of integers: where to write
                      max_len = 80 # max. length of line
                      ):
    """Constructs string containing 'words' starting at 'postitions'"""
    if type(words) != list:
        words = [words]
        positions = [positions]
    if len(words) != len(positions):
        raise ValueError("'words' and 'positions' must be the same length.")
    s = ''
    for k in range(len(words)):
        s = ' '*(positions[k] - len(s)) + str(words[k])
    if min(positions[-1] + len(str(words[-1])), max_len) < len(s):
        raise ValueError("Too long words.")
    f.write(s + '\n')

