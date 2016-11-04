# Built-in modules #
import os, re, sys, warnings
#
# Third party modules #
import numpy as np
import sh, shutil
from subprocess import call
import matplotlib.pyplot as plt
#
# Short forms #
npcat = np.concatenate
#

class Main(object):
    """
    'God-ish class' which can be used as a user interface. The initialization
    takes an ENDF file, sets up resonance file (self.resonance_file), creates 
    a SammyRunner instance (self.sr) and a PseudoExperiment instance 
    (self.experiment) using self.resonance_file. self.sr or
    self.experiment may be modified (to change settings). Then, 
    self.fit(out_paths) can be called to finalize the retroactive method by 
    fitting parameters to self.experiment using self.resonance_file as 
    starting guess. The resulting ResonanceFile instance is stored as 
    self.resonance_out and the resulting files are stored at the entries of 
    out_paths (list of length 2).
    
    Example:
    >> import retroactipy as ret
    >> main = ret.Main('Ni059.endfb71')
    >> main.fit(['out.mf2', 'out.mf32'])
    >> self.resonance_out.path
       '/home/petter/pythonWork/course/retroactipy/out.mf2'
    >> self.resonance_out.cov.path
       '/home/petter/pythonWork/course/retroactipy/out.mf32'
    >>
    
    The resulting ENDF formatted resonance and resonance covariance files are 
    stored at 'out.mf2' and 'out.mf32', respectively in this example.
    """
    
    def __init__(self, endf_path):
        """Constructor for Main objects"""
        self.endf_file = CompleteEndfFile(endf_path)
        self.resonance_file = self.endf_file.create_resonances(
            '_retroactipy.mf2', '_retroactipy.mf32')
        self.sr = SammyRunner()
        self.add_experiment()
        
        
    def add_experiment(self,
                       npoints = 5000,
                       nsections = 10,
                       least_E_min = 1e-2,
                       ):
        E_min = max(self.resonance_file.min_energy, least_E_min)
        E_max = self.resonance_file.max_energy
        E_section_lims = np.logspace(np.log10(E_min), np.log10(E_max), 
                                     nsections + 1)
        time_of_flight = npcat([np.linspace(1/np.sqrt(E_section_lims[i]),
            1/np.sqrt(E_section_lims[i+1]), npoints/nsections + 1)[:-1]
            for i in range(len(E_section_lims)-1)])
        energies = 1/time_of_flight**2
        self.experiment = PseudoExperiment(
            self.sr.reconstruct(self.resonance_file, energies))
        
    def fit(self, out_paths):
        self.resonance_out = self.sr.fit(self.resonance_file, 
            self.experiment, out_paths)
        return self.resonance_out


class EndfFile(object):
    """
    Class describing EndfFile objects. Only contains a path, which is error
    checked at initialization. The class could be expanded with attributes and
    methods that are common for all ENDF formatted files.
    """
    
    def __init__(self, path):
        """Constructor with error checks for the path"""
        path = absolute_path(path)
        assert path != ''
        if not os.path.exists(path):
            raise ValueError('Path %s does not exist' % path)
        assert not os.path.isdir(path)
        self.path = path

        
class CompleteEndfFile(EndfFile):
    """Class describing the subset of ENDF files which are complete."""
    
    def __init__(self, path):
        EndfFile.__init__(self,path)
    
    def create_resonances(self, path, pathcov): 
        """Creates (and returns) ResonanceFile instance and produces the 
        corresponding file(s) containing the subsection of the 
        ENDF file with resonance parameters and (if available) 
        covariances"""
        try: 
            # Check if parsing for nuclide is already done #
            self.nuclide
            #
        except AttributeError:
            # Add nuclide attribute #
            self.nuclide = self.create_nuclide()
            #
        # Open files for reading and writing #
        f = open(self.path,'r')
        fres = open(path,'w')
        frescov = open(pathcov,'w')
        #
        # Booleans indicating if certain sections are found #
        found_res = False
        found_cov = False
        #
        while True: # To enable additional reading in loop
            line = f.readline()
            if line == '': break # EOF
            # Generate EndfLineParser for easy access #
            el = EndfLineParser(line)
            #
            # Look for MF = 2 (resonance section) #
            if el['mf'] == 2 and el['mt'] == 151:
                found_res = True
                # Copy resonance section
                while el['mf'] == 2 and el['mt'] == 151:
                    fres.write(el.line)
                    el = EndfLineParser(f.readline())
                fres.write(el.line) # one more line shall be included
                #
                # Create corresp. ResonanceFile object #
                fres.flush() # initialization needs reading
                res_file = ResonanceFile(path, self.nuclide)
                #
            #
            # Look for MF = 32 (covariance section) #
            elif el['mf'] == 32 and el['mt'] == 151:
                found_cov = True
                # Copy covariance section #
                while el['mf'] == 32 and el['mt'] == 151:
                    fres.write(el.line)
                    el = EndfLineParser(f.readline())
                fres.write(el.line) # one more line shall be included
                #
                # Create corresp. ResonanceCovFile object #
                res_file.cov = ResonanceCovFile(pathcov)
                #
            elif found_res and found_cov:
                # Both are found; no need to continue #
                break
                #
            #
        # Close files
        for fi in [f, fres, frescov]: fi.close()
        #
        # If no MF = 2 found, raise exception
        if not found_res:
            os.remove(path)
            raise ValueError('No resonance "file" in %s' % self.path)
        #
        # If covariance file empty: remove
        if not found_cov: os.remove(pathcov)
        #
        
        return res_file
            
                        
    def create_nuclide(self):
        """Parses head of the file and calls create resonances,
        returns Nuclide instance"""
        with open(self.path) as f:
            f.readline() # Skip zeroth line
            # Line 1 #
            el = EndfLineParser(f.readline())             
            za = int(el[0]) # float combining atomic number and mass
            neutron_mass = 1.00866491588 # Natural constant
            weight = el[1] * neutron_mass
            atomic_number = za / 1000
            mass = za - atomic_number * 1000
            mat_nbr = el['mat']
            #
            # Line 2
            el = EndfLineParser(f.readline())
            state = el[2]
            #
            # Go to 5th line
            for k in range(3,6): el = EndfLineParser(f.readline()) 
            #
            element = el[0].split('-')[1]
            # Initialize nuclide #
            self.nuclide = Nuclide(atomic_number, element, mass, state)
            #
            # Add more information for easy access #
            self.nuclide.mat_nbr = mat_nbr
            self.nuclide.weight = weight
            self.nuclide.channel_radius = 0.
            #
            return self.nuclide
     
        
class ResonanceFile(EndfFile):
    """
    Class describing ENDF formatted files which contain the resonance part 
    only (MF = 2, MT = 151), the file can be created by method in the
    CompleteEndfFile class.
    """
    
    def __init__(self, path, nuclide):
        EndfFile.__init__(self, path)
        self.nuclide = nuclide
        # Parse to get basic information #
        with open(self.path,'r') as f:
            f.readline()
            f.readline()
            el = EndfLineParser(f.readline())
            self.min_energy, self.max_energy = el[:2]
        #
            
    def plot(self, xs_type = 'total', **kwargs):
        """
        Reconstructs the cross sections (if necessary) and plots the cross 
        section of type 'xs_type', which can be either mt number or reaction
        type as defined in SAMMY.
        """
        # Check if already reconstructed, store as local PseudoExperiment #
        # instance for easy access.
        try: xss = PseudoExperiment(self.cross_sections)
        #
        except AttributeError:
            # Reconstruct if necessary
            sr = SammyRunner()
            energies = np.logspace(np.log10(max(self.min_energy, 1e-2)), 
                np.log10(self.max_energy), 2000)
            self.cross_sections = sr.reconstruct(self, energies)
            xss = PseudoExperiment(self.cross_sections)
            #
        # Plot using method of CrossSection class
        xss[xs_type].plot(**kwargs)
            
        
class ResonanceCovFile(EndfFile):
    """
    Class describing ENDF formatted files which contain the resonance 
    covariance part only (MF = 32, MT = 151), the file can be created by 
    method in thev CompleteEndfFile class. More methods may be added.
    """
    def __init__(self, path):
        EndfFile.__init__(self, path)



class Nuclide(object):
    """Class describing the basic features of a nuclide"""
    
    def __init__(self, atomic_number, element, mass, state):
        """Constructor of Nuclide objects"""
        self.atomic_number = int(atomic_number)
        self.element = element
        self.mass = int(mass)
        self.state = int(state)
        

class CrossSection(object):
    """Class describing a cross section."""
    
    def __init__(self, 
                 mt, # Identifying so called mt number (ENDF jargon)
                 energies, # Neutron energy grid (numpy array)
                 cross_section): # Cross section at 'energies' (numpy array)
        """Constructor of CrossSection objects"""
        self.mt = int(mt) 
        self.energies = np.array(energies)
        self.cross_section = np.array(cross_section)
        self.experiment = None
        # Translate mt number into other forms for easy access #
        mts = np.array([1, 2, 102, 18]) # Possible mt numbers
        xs_types = ['total', 'elastic', 'capture', 'fission']
        xs_strings = ['(n,tot)', '(n,el)', '(n,$\gamma$)', '(n,f)']
        idx = np.where(mt == mts)[0]
        self.sammy_type = xs_types[idx]
        self.string = xs_strings[idx]
        #
        
    def __str__(self):
        """Returns pretty string describing the cross section type"""
        return self.string
        
    def plot(self, **kwargs):
        """Plots cross section as function of energy"""
        plt.loglog(self.energies, self.cross_section, **kwargs)
        ax = plt.gca()
        ax.set_xlabel('Energy [eV]')
        ax.set_ylabel('%s cross section [b]' % str(self))
        ax.grid(True)
        
class PseudoExperiment(object):
    """
    Class describing a set of pseudo-experimental cross section data. This
    includes CrossSection objects including uncertainty information, and can 
    include some more information about the imagined experiment.
    """
    
    def __init__(self, 
                 cross_sections, # List of CrossSection objects
                 thickness = 1e-3 # Assumed sample thickness [atoms/b = b^(-1)]
                 ):
        """Constructor for PseudoExperiment objects"""
        self.cross_sections = cross_sections
        # Let each CrossSection instance know about its owner #
        for xs in cross_sections:
            xs.experiment = self
        #
        self.thickness = thickness
        # Generate uncertainty information #
        self.generate_unc()
        #
        
    def __getitem__(self, idx):
        """
        Uses mt number or type (at least three first chars) for indexing.
        A list of such indices is accepted, but not slices.
        """
        tp = type(idx)
        if tp == list:
            # Return list corresponding to lis of indices #
            answer = []
            for i in idx:
                answer.append(self[i])
            return answer
            #
        elif tp == int:
            # Look for CrossSection with this mt number #
            for xs in self.cross_sections:
                if xs.mt == idx:
                    return xs
            #
            raise IndexError(
                'mt = %s not found in PsedoExperiment instance.' % str(idx))
        elif tp == str:
            # Look for CrossSection with this reaction type #
            for xs in self.cross_sections:
                if xs.sammy_type[:max(3,len(idx))] == idx:
                    return xs
            #
            raise IndexError(
                'xs type = %s not found in PsedoExperiment instance.' % idx)
        else:
            # If idx not of any of the above types:
            raise ValueError('%s type not allowed for indexing.' % str(tp))
            #
        
        
    def generate_unc(self):
        """
        Function which generates assumed uncertainty information for the   
        experiment. Contains hard-coded assumptions, which should be avoided.
        Should add special treatment of transmission experiments.
        """
        for xs in self.cross_sections:
            # 10 % statistical unc., will have little impact #
            xs.statistical_unc = xs.cross_section*.10 
            #
            # Create dictionary with uncertainty information for SAMMY. #
            # Energy resolution information can be added here. #
            xs.unc = {}
            xs.unc['thickness'] = self.thickness*0.02
            xs.unc['normalization'] = 0.005 # relative
            # List with background unc.s: constant backg unc as first element #
            # More background types can be added. #
            xs.unc['backgrounds'] = [0.01] # [barn]
            #
            xs.unc['detector efficiency'] = 0.01 # relative
            #
        
    def plot(self, xs_type = 1, **kwargs):
        """
        Error bar plot of the experiment, the uncertainties are statistical 
        only.
        """
        xs = self[xs_type]
        lowlim = 1e-8
        plt.errorbar(xs.energies, xs.cross_section,yerr = np.array([np.minimum(
            xs.statistical_unc,xs.cross_section-lowlim),xs.statistical_unc]),
            **kwargs)
        ax = plt.gca()
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_xlabel('Energy [eV]')
        ax.set_ylabel('%s cross section [b]' % str(xs))
        ax.grid(True)


class EndfLineParser(object):
    """
    Class which is initialized using a line from an ENDF file, and then can be
    used to call for certain parts of that line, using square brackets. 
    Indices of normally formatted lines (or oddly formatted using self.breaks)
    can be used, or certain keywords. Lists combining these two, and slices, are
    also allowed. 
    
    Possible keywords:
    'mf' or 'file' for mf number
    'mt' or 'tape' for mt number
    Anything starting with 'mat' for material number
    Anything starting with 'num' or 'nbr' for 'record number' (ENDF jargon)
    """
    
    def __init__(self, 
                 line, # Line to parse
                 breaks = None # Columns breaks (if None, set to ENDF default)
                 ):
        self.line = line
        self.default_breaks = npcat((np.array(range(0,7,1))*11,[70,72,75,80]))
        # Custom column breaks are rare; assume default
        if breaks == None:
            self.breaks = self.default_breaks
        #
        
    def __getitem__(self,idx):
        """
        Method to return elements of the EndfLineParser instance. Accepts 
        indices (in range), certain keywords, lists of such, or slices.
        """
        tp = type(idx)
        if tp == int:
            # Always ends here (possibly after recursion)
            # Returns integer if nothing indicates float. If float,
            # ENDFs weird format for floating point numbers needs to be treated
            # If impossible to parse into integer or float, return str (can
            # be allowed).
            idx = idx % len(self) # to get correct idx + 1
            # Get the str
            word = self.line[self.breaks[idx]:self.breaks[idx+1]] 
            #
            if word.strip().isdigit():
                # Integer if "possible"
                return int(word)
                #
            try:
                # Try to return float #
                # Identify +/- signs that are powers: #
                m = re.findall('([0-9.]+)([+-][0-9]+)',word)
                if m != []: 
                    # power in ENDF format found
                    # reconstruct with E for power:
                    return float(m[0][0] + 'E' + m[0][1])
                    #
                #
                elif re.search('[^ \t\n\r\f\v]',word) == None:
                    # Empty entry => float(0)
                    return 0.
                else:
                    return float(word)
                #
            except ValueError:
                # if not a number, return string
                return word
            #
        elif tp == list:
            # Return list with corresponding entries #
            answer = []
            for i in idx:
                answer.append(self[i])
            return answer
            #
        elif tp == slice:
            # Set up corresponding list, recurse with this list #
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
            #
        elif tp == str:
            # String should be one of the allowed keywords #
            self.breaks = self.default_breaks
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
            #
        else:
            raise TypeError('%s invalid as index.' % str(tp))
        
    def __len__(self):
        """
        Defines the length as the number of entries in the parsed line.
        """
        return len(self.breaks)-1
        

class SammyRunner(object):
    """
    Class collecting methods and settings (attributes) for running SAMMY.
    Structure-wise uninteresting, but handles the most important methods.
    """
    
    def __init__(self, 
            executable = 'sammy', 
            # SAMMY settings #
            data_format = 'twenty', # Format used in reading writing data
            temperature = 0., # This should be changed in future work
            default_thickness = 0.003, # Sample thickness if not in experint
            baby_steps = True, 
            lcomp = 2,
            ndigits_f32 = 6,
            include_correlations = True,
            #
            cleanup = True, # Should temporary files be removed?
            # Convergence criteria in GLSQ wrapper #
            convergence_dict = dict(
                        max_iter = 200, # Max nbr of iterations
                        length = 10, # Number of chi-sq to consider
                        tol = 0.01, # Reduced chi-sq tol
                   ),
            #
            ):
    
        # Settings are stored as attributes #
        self.executable = executable
        if not sh.which(self.executable):
            warnings.warn(
                "Sammy executable '%s' not found.\n" % self.executable + 
                "Please install SAMMY or change executable.")
        self.data_format = data_format
        self.temperature = temperature
        self.default_thickness = default_thickness
        self.baby_steps = baby_steps
        self.lcomp = lcomp
        self.ndigits_f32 = ndigits_f32
        self.include_correlations = include_correlations
        self.default_xs_type = 'total'
        self.cleanup = cleanup
        self.convergence_dict = convergence_dict
        self.convergence_dict['min_iter'] = self.convergence_dict['length']
        #
        
    def reconstruct(self, 
    				resonance_file, 
		    		energies # np.array with energies, or list of such arrays
		    		):
        """
        Reconstructs the cross sections using parameters in resonance_file 
        object at the provided energies, and using the 'settings' of the 
        SammyRunner object. As of now, the reconstructions are done at SAMMYs 
        default energies and interpolated, this should be changed and should 
        enable Doppler Broadening (self.temperature > 0.). Fission shall be 
        added.
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
        tempinp = temp_file_gen('Sammy_reconstruct','inp')
        temppar = temp_file_gen('Sammy_reconstruct','par')
        templst = temp_file_gen('Sammy_reconstruct','lst')
        #
        # Generate input files from ENDF
        self.endf2inp_par_ndf(resonance_file, [tempinp, temppar])
        #
        # Prepare actual reconstruction run by modifying input #
        newkeylines = ['do not solve bayes equations',
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
        
    def fit(self, resonance_file, experiment, out_paths):
        """
        Produces SAMMY input and runs SAMMY 'without prior',
        returns a new ResonanceFile instance with resonance parameters fit 
        to the pseudo-experiment. The fit is so far simplified such that 
        'true GLSQ' is used for (n,tot) and then Bayesian updating is used 
        for (n,gamma). This would be equivalent to GLSQ all the way   
        if there were no non-linearities and no correlations between the
        different reactions. All this has to be developed further later 
        (outside course).
        """
        # Set up temporary file names #
        inp = temp_file_gen('Sammy_fit','inp')
        par = temp_file_gen('Sammy_fit','par')
        cov = temp_file_gen('Sammy_fit','cov')
        ndf = temp_file_gen('Sammy_fit','ndf')
        parout = temp_file_gen('Sammy_fit','out.par')
        covout = temp_file_gen('Sammy_fit','out.cov')
        #
        # Construct SAMMY input using resonance_file and information about the #
        # 'experiment' #
        self.endf2inp_par_ndf(resonance_file, [inp, par, ndf], 
            experiment[1], flag_all = True)
        #
        # Change from MLBW formalism if this was in original file. #
        # Reich-Moore will be used instead, which is recommended. #
        self.modify_inp(inp, keyremove = ['mlbw formalism is wanted'])
        #
        # Fit to total cross section data without prior #
        message = self.g_least_squares(inp, par, experiment['total'],
            parout, covout)
        shutil.move(parout, par)
        shutil.move(covout, cov)
        #
        # Check if convergence was reached. Otherwise, something is bad. #
        if message[:len('Did not converge')] == 'Did not converge':
            raise RuntimeError(message)
        #
        # Perform a Beyesian update using capture data
        self.bayesian([inp, par, cov], experiment['capture'], [parout, covout])
        #
        # Construct ENDF formatted files from output #
        self.inp_par_ndf_cov2endfs([inp, parout, ndf, covout], out_paths)
        #
        # Include ENDF file paths in ResonanceFile instance to return
        resonance_file_out = ResonanceFile(out_paths[0], resonance_file.nuclide)
        resonance_file_out.cov = ResonanceCovFile(out_paths[1])
        #
        # Clean up
        if self.cleanup:
            for p in [inp, par, cov, ndf, parout, covout]: os.remove(p)
        #
        return resonance_file_out
        
    def endf2inp_par_ndf(self, 
                resonance_file, # Implicitly containing resonance parameters etc
                out_paths, # List of desired paths for output: .inp, .par, .ndf
                xs_measurement = None, # Possible PseudoExperiment instance
                flag_all = False, # If all res. parameters should be flagged
                ):
        """
        Given a ResonanceFile instance and a list with the paths to desired .inp
        .par and .ndf (optional) files, the function uses SAMMY to produce 
        these files. These files can then be ued as input to other SAMMY runs.
        """
        nucl = resonance_file.nuclide # For easy access
        # Set name of temporary file #
        tempinp = temp_file_gen('Sammy_endf2par', 'inp')
        # 
        # Set up specific input lines to SAMMY #
        keylines = [
                    'do not solve bayes equations',
                    'input is ENDF/B file MAT=%i' % nucl.mat_nbr,
                    'use energy range from ENDF/B file 2',
        ]
        #
        # Add lines if .ndf file is wanted #
        if len(out_paths) > 2: 
            keylines.insert(0,'automatic ndf file creation')
            keylines.append('ndf file is in key word format')
        #
        # Generate title for SAMMY input of a max. length of 80 chars #
        title = 'SAMMY inp (generated from) for ENDF '
        if len(title) + len(resonance_file.path) > 80: 
            pathstr = '...' + resonance_file.path[-(77 - len(title)):]
        else: pathstr = resonance_file.path
        title += pathstr
        #
        # Construct the input file #
        self.construct_inp(tempinp, nucl, keylines, title, xs_measurement)
        # 
        # Set up which output from SAMMY is desired #
        desired = ['SAMNDF.INP', 'SAMNDF.PAR', 'SAMNDF.NDF'][:len(out_paths)]
        #
        # Run SAMMY with tempinp and resonance_file.path as input files
        self.run([tempinp, resonance_file.path, 'dummy.dat'],
            desired, out_paths)
        #
        # Modify obtained input for use in other steps
        keyremove = ['do not solve bayes equations',
                     'automatic ndf file creation']
        keylines = []
        if flag_all: keylines.append('flag all resonance parameters')
        try: thickness = xs_measurement.experiment.thickness
        except AttributeError: thickness = None
        self.modify_inp(out_paths[0], keylines, keyremove = keyremove, 
            new_thickness = thickness)
        #
        if len(out_paths) > 2:
            # If .ndf file is produced, modify its content:
            # Change LCOMP to self.lcomp and set NDIGITS to self.ndigits_f32
            tempndf = temp_file_gen('Sammy_endf2par', 'ndf') 
            if len(out_paths) > 2:
                with open(out_paths[2],'r') as f1:
                    with open(tempndf,'w') as f2:
                        for line in f1:
                            if line[:5].lower() == 'lcomp':
                                f2.write('Lcomp = %i\n' % self.lcomp)
                                if self.lcomp == 2:
                                    f2.write('Ndigits = %i\n' 
                                        % self.ndigits_f32)
                            else:
                                f2.write(line)
        
            shutil.move(tempndf, out_paths[2])
        #
        # Clean up #
        if self.cleanup:
            os.remove(tempinp)
        #
        
        
    def inp_par_ndf_cov2endfs(self, input_paths, out_paths):
        """
        Given a list with SAMMY .inp, .par, .ndf, and .cov files, this method
        produces a resonance file and a resonance covariance file of  ENDF 
        format at out_paths.
        
        Note for later work (outside course):
        File 32 'alone is not capable...' See SAMMY-8 manual rev 7 p. 532.
        """
        # Specific input to SAMMY #
        keylines = [
                    'do not solve bayes equations',
                    'ENDF/B-VI file 2 is wanted',
                    'ndf file is in key word format',
                    'put covariance matrix into endf file 32'
                   ]
        #
        # Set up and modify temporary input file #
        tempinp = temp_file_gen('Sammy_inps2endfs','inp')
        shutil.copy(input_paths[0], tempinp)
        input_paths[0] = tempinp
        self.modify_inp(input_paths[0], keylines)
        #
        # Run SAMMY
        self.run(input_paths, ['SAMMY.NDF', 'SAMMY.N32'], out_paths)
        #
        # Clean up
        if self.cleanup: os.remove(tempinp)
        #
                
    def g_least_squares(self, 
             sammy_inp, sammy_par, # assumed out_paths of endf2inp_par_ndf()
             xs_measurement, # PseudoExperiment instance
             out_par, out_cov
             ):
        """
        Performs GLSQ fit to one type of cross section using SAMMY itertively 
        until convergence. 
        """
        # Set up out paths #
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
        try: # Combined with finally to ensure moving back to original dir #
            # Set file names #
            dat_path =  'data.dat'
            inp_path = 'inp.inp'
            #
            # Construct .dat file #
            data = np.array([xs_measurement.energies, 
                            xs_measurement.cross_section,
                            xs_measurement.statistical_unc]).T
            self.write_data(data, dat_path)
            #
            # Modify input #
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
            if self.include_correlations:
                # Add experimental uncertainties to inp file #
                self.exp_unc_into_inp(inp_path, xs_measurement)
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
                # Set input to self.run() #
                output_names = [next_par, output, cov]
                inp_paths[1] = current_par
                if k == 0:
                    shutil.copy(originals[1], current_par)
                #
                # Run SAMMY #
                self.run(inp_paths, desired_outputs, output_names)
                #
                # Find number of points and chi-squared in output #
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
            else: # Reached if no break statement
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
        # Set up temporary files #
        tempinp = temp_file_gen('Sammy_bayesian','inp')
        tempdata = temp_file_gen('Sammy_bayesian','dat')
        shutil.copy(sammy_inputs[0], tempinp)
        sammy_inputs[0] = tempinp
        # 
        # Construct .dat file from xs_measurement #
        data = np.array([xs_measurement.energies, 
                        xs_measurement.cross_section,
                        xs_measurement.statistical_unc]).T
        self.write_data(data, tempdata)
        #
        # Include .dat file in list of SAMMY inputs #
        sammy_inputs.insert(2,tempdata)
        #
        # Change reaction type in input file #
        self.modify_inp(sammy_inputs[0], new_type = xs_measurement.sammy_type)
        #
        if self.include_correlations:
            # Add experimental uncertainties to inp file #
            self.exp_unc_into_inp(sammy_inputs[0], xs_measurement)
            #
        # Run SAMMY #
        self.run(sammy_inputs, ['SAMMY.PAR','SAMMY.COV'], 
                 output_names = sammy_outputs)
        #
        # Clean up #
        if self.cleanup:
            for p in [tempinp, tempdata]: os.remove(p)
        #
        
            
    def construct_inp(self, path, nuclide, keylines, title = None, 
            xs_measurement = None):
        """
        Given certain information, in 'self', 'nuclide', and 'keylines',
        and possibly also within 'xs_measurement', the method constructs a 
        basic .inp file at 'path'.
        """
        # Add data format to keylines - should always be there #
        if not self.data_format in keylines:
            keylines = [self.data_format] + keylines 
        #
        # Get thickness and reaction type from xs_measurement if possible, #
        # otherwise use defaults #
        try: thickness = xs_measurement.experiment.thickness 
        except: thickness = self.default_thickness
        try: xs_type = xs_measurement.sammy_type
        except: xs_type = self.default_xs_type
        #
        # Set title #
        if title == None:
            title = 'SAMMY inp file for %s-%s' % (nuclide.element, nuclide.mass)
        #
        with open(path,'w') as f:
            # Write to the .inp file #
            f.write(title[:80] + '\n')
            write_to_position(f, ['%s%i' % (nuclide.element, nuclide.mass), 
                nuclide.weight], [0, 10])
            for keyline in keylines + ['']:
                f.write(keyline + '\n')
            write_to_position(f, self.temperature, 0)
            write_to_position(f, [nuclide.channel_radius, thickness],[0, 10])
            write_to_position(f, xs_type, 0)
            #

    def modify_inp(self, path, keylines = [], keyremove = [], 
                         new_type = None, new_thickness = None):
        """
        Method to modify an existing .inp file. Can add and remove keylines,
        change reaction type and change sample thickness.
        """
        # Set up temporary .inp file name #
        tempnew = temp_file_gen('Sammy_modify_inp', 'inp')
        #
        # Counter for line number after the keylines are finished #
        linenbr_after_keylines = 0
        #
        with open(path,'r') as f1:
            with open(tempnew,'w') as f2:
                for line in f1:
                    if line.strip() == '':
                        # Found end of current keylines, add keylines here #
                        for keyline in keylines:
                            f2.write(keyline + '\n')
                        #
                        # Start counting lines after keylines #
                        linenbr_after_keylines = 1
                        #
                    if (new_type != None and line[:3].lower() in 
                               ['cap','tot','ela','tra'] and 
                               linenbr_after_keylines == 4):
                        # Line with reaction found and shall be replaced #
                        f2.write(new_type + '\n')
                        #
                    elif new_thickness != None and linenbr_after_keylines == 3:
                        # Line with thickness found and shall be replaced #
                        lstline = list(line)
                        lstline[10:20] = list(str(new_thickness).rjust(10))
                        f2.write(''.join(lstline))
                        #
                    else:
                        # Copy line if shall not be removed #
                        if not any([line.lower()[:len(keyrem)] == keyrem.lower()
                                for keyrem in keyremove + keylines]):
                            f2.write(line)
                        #
                    if linenbr_after_keylines > 0:
                        linenbr_after_keylines += 1
        # Move the modified file to original path
        shutil.move(tempnew, path)
        #
        
    def exp_unc_into_inp(self, inp, xs_measurement):
        """
        Method which writes experimental systematic uncertainty information
        into .inp file. Handles uncertainty in sample thickness, background of 
        SAMMYs various types, normalization, and a simple case of detector 
        efficiency.
        """
        # Lines in .inp used to announce data reduction parameters #
        # and their uncertainties #
        announce_lines = [
                          'broadening parameters may be varied',
                          'normalization and background are next',
                          'detector efficiencies'
        ]
        #
        # Length of these lines that are necessary to include #
        annonuce_necessary_lengths = [5, 5, 5]
        #
        # Set up temporary file name #
        tempinp = temp_file_gen('Sammy_exp_unc_into_inp','inp')
        #
        with open(tempinp, 'w') as f2:
            with open(inp, 'r') as f1:
                # Copy from original input until end or until data reduction #
                # parameters are reached (which must be at the end) # 
                for line in f1:
                    if any([line[:annonuce_necessary_lengths[i]].lower() == 
                            announce_lines[i][:annonuce_necessary_lengths[i]] 
                            for i in range(len(announce_lines))]):
                        break
                    f2.write(line)
                #
            # Check if list of uncertainties to ignore is added to SammyRunner #
            # instance (for debugging) #
            try: skip_list = self.skip_list
            #
            # Otherwise, start with empty skip_list #
            except AttributeError: skip_list = []
            #
            xsmunc = xs_measurement.unc # Short form
            # Loop over the unc dictionary with systematic uncertainties
            for key, content in xsmunc.items():
                if key in skip_list:
                    pass
                elif key == 'thickness':
                    # Thickness uncertainty #
                    f2.write(announce_lines[0] + '\n')
                    write_to_position(f2, 
                        [xs_measurement.experiment.thickness, 3], [20, 65])
                    write_to_position(f2, content, 20)
                    f2.write('\n')
                    #
                elif key == 'normalization': 
                    # Normalization uncertainty (and background if available)
                    f2.write(announce_lines[1] + '\n')
                    if not 'backgrounds' in xsmunc.keys():
                        xsmunc['backgrounds'] = []
                    values = [1.] + [0.]*len(xsmunc['backgrounds'])
                    value_pos = range(0,60,10)[:len(values)]
                    flags = [3]*len(values)
                    flags[0] = 3*(content > 0)
                    flag_pos = range(60,72,2)[:len(flags)]
                    write_to_position(f2, values + flags, 
                        value_pos + flag_pos)    
                    unc = [content] + xsmunc['backgrounds']
                    write_to_position(f2, unc, value_pos)
                    f2.write('\n')
                    skip_list.append('backgrounds') # Ignore 2nd time
                    #
                elif key == 'backgrounds':
                    # Is included together with normalization
                    if not 'normalization' in xsmunc.keys():
                        xsmunc['normalization'] = 0.
                    #
                elif key == 'detector efficiency':
                    f2.write(announce_lines[2] + '\n')
                    write_to_position(f2, [1., content, 3], [0, 10, 20])
                    f2.write('\n')
                else:
                    raise ValueError('Uncertainty key %s not implemented' 
                                      % str(key))        
                    # More uncertainty types need to be implemented
        # Move modified .inp to path of original
        shutil.move(tempinp, inp)
        #
        
    def run(self, input_paths, desired_outputs, output_names = None):
        """
        Given a list of paths to input files (possibly dummies), desired output
        files and desired names of these output files, the method constructs
        a shell script which runs SAMMY and answers the prompt such that it is
        fed with the input, and then obtains the desired output files at their 
        respective paths.
        """
        # Set up and move to temporary directory #
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
        #
        try: # Combined w. finally to ensure moving back to original directory #
            # Construct shell script 'run_file' #
            run_file = 'run_file'
            with open(run_file, 'w') as f:
                f.write('#! /bin/sh\n\nsammy<<EOF\n')
                for path in input_paths + ['','EOF']:
                    f.write(path + '\n')
            #
            # Change permissions #
            sh.chmod('755', run_file)
            #
            # Run shell script running SAMMY and answering prompts #
            call('./%s > terminal.out 2>&1' % run_file, shell = True)
            #
            # Rename the desired outputs #
            for k in range(len(desired_outputs)):
                shutil.move(desired_outputs[k],output_names[k])
            #
        # Move back to original directory even if error has been raised #
        finally:
            os.chdir(origdir)
        #
        # Clean up #
        if self.cleanup: sh.rm('-rf',tempdir)
        #
            
    def read_data(self, path):
        """
        Reads data in 'path' using format in 'self.data_format', returns 
        numpy array. Only SAMMYs 'format twenty' is allowed, but more formats
        can easily be added.
        """
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
    
def temp_file_gen(base, # Base name
                  extension = '', # File extension
                  directory = False, # Make directory instead of empty file
                  always_first = 'temp_'
                  ):
    """
    Generates name of temporary file including a number and a possible file
    extension, creates empty file.
    """
    base = re.sub('[^/]*$','',base) + always_first + re.sub(
                                                    '.*/','',base).strip('_')
    extension = '.'[:(len(extension)*(extension[:1] != '.'))] + extension
    k = 0
    while True:
        attempt = base + '_' + str(k).zfill(3) + extension
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
        s += ' '*(positions[k] - len(s)) + str(words[k])
    if min(positions[-1] + len(str(words[-1])), max_len) < len(s):
        raise ValueError("Too long words.")
    f.write(s + '\n')

