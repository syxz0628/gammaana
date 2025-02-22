import numpy as np
import argparse
import pymedphys
import pydicom
from datetime import datetime, date, time
import time
import os
import nrrd
import sys
import matplotlib.pyplot as plt
import sys
class class_gammaanalysis():
    def __init__(self):
        self.header=[]
    def fun_readDoseNRRD(self,filepath, **kwargs):
        directory = filepath[:filepath.rfind("/") + 1]
        temppath = directory+"MakeContour_tempFileWoSpaceDir.nrrd"
        space_directions = []
        space_origins = []
        dose_data = []
        dose_data_size = []

        with open(filepath) as inFile:
            with open(temppath, "w") as outFile:
                for line in inFile:
                    if line.find("space directions") != -1:
                        try:
                            v = line.split()
                            space_directions = [float(v[2][1:-1]), float(v[6][:-1]), float(v[-1][:-1])]
                            continue
                        except:
                            pass
                    if 'space origin' in line:
                        v = line.split()
                        s = v[2].split(',')
                        space_origins = [int(s[0][1:]), int(s[1]), int(s[2][:-1])]
                    if line.find("data file") != -1:
                        v = line.split(' ')
                        new = directory + v[-1]
                        line = v[0] + ' ' + v[1] + ' ' + new
                    outFile.write(line)

        data, self.header = nrrd.read(temppath)
        os.system("rm " + temppath)
        data = np.array(data)

        # do = plt.imshow(np.rot90(dose_masked[slice_x_dose, slice_y_dose, slice_z_dose]), origin='lower',
        #                 extent=[0., dose_masked[slice_x_dose, slice_y_dose, slice_z_dose].shape[0] * pixel_size_dose[0], 0.,
        #                         dose_masked[slice_x_dose, slice_y_dose, slice_z_dose].shape[1] * pixel_size_dose[1]],
        #                 cmap=cmap, alpha=dose_alpha, norm=norm)
        data = np.resize(data, (self.header['sizes'][0], self.header['sizes'][1], self.header['sizes'][2]))
        if not np.any(space_directions):
            try:
                space_directions = [self.header['space directions'][0, 0], self.header['space directions'][1, 1],
                                    self.header['space directions'][2, 2]]
            except:
                print('Could not read the space directions from the nrrd file')
                raise ValueError

        if kwargs.get('overwrite', False) and np.any(dose_data):  ## do not override if there is nothing to override
            dose_file = kwargs.get('dose_file', -1)
            dose_data[dose_file] = data
            dose_data_size[dose_file] = space_directions
        else:
            #print(f'The dose cube has dimensions {data.shape}')
            dose_data.append(data)
            dose_data_size.append(space_directions)

        axes_array = [[], [], []]
        j = space_origins[0]
        for i in range(0, self.header['sizes'][0]):
            j = j + space_directions[0]
            axes_array[0].append(float(j))

        j = space_origins[1]
        for i in range(0, self.header['sizes'][1]):
            j = j + space_directions[1]
            axes_array[1].append(float(j))

        j = space_origins[2]
        for i in range(0, self.header['sizes'][2]):
            j = j + space_directions[2]
            axes_array[2].append(float(j))

        axes_array_np = tuple(axes_array)
        return axes_array_np, data

    def fun_gamma_analysis(self,additionalinfo,dose1, dose2, dosediscrit, cuoff, maxdose, interfra, maxgamma, fraction, saveresultas,pronecase,moreinfo):
        if dose1[dose1.rfind('.'):] == '.dcm' and dose2[dose2.rfind('.'):] == '.dcm':
            reference = pydicom.dcmread(dose1)
            evaluation = pydicom.dcmread(dose2)
            axes_reference, dose_reference = pymedphys.dicom.zyx_and_dose_from_dataset(reference)
            axes_evaluation, dose_evaluation = pymedphys.dicom.zyx_and_dose_from_dataset(evaluation)
        elif dose1[dose1.rfind('.'):] == '.nrrd' and dose2[dose2.rfind('.'):] == '.nrrd':
            axes_reference, dose_reference = self.fun_readDoseNRRD(dose1)
            axes_evaluation, dose_evaluation = self.fun_readDoseNRRD(dose2)
        else:
            print("wrong dose file (both ref and com should be nrrd or dcm dose files only), check input.")
            sys.exit()

        dose_evaluation_fx=dose_evaluation*float(fraction)
        if pronecase:
            dose_ref_temp = np.fliplr(dose_reference[:, :, ])
            dose_reference = np.flipud(dose_ref_temp[:, :, ])
            write2file=dose2[:dose2.rfind('.')]+'_rewrite_ref.nrrd'
            nrrd.write(write2file, dose_reference, self.header)
            write2file=dose2[:dose2.rfind('.')]+'_rewrite_com.nrrd'
            nrrd.write(write2file, dose_evaluation_fx, self.header)

        # p = dose_reference[:, :, 70]
        # plt.imshow(np.fliplr(p), origin='lower')
        # plt.show()
        # p = dose_evaluation_fx[:, :, 70]
        # plt.imshow(np.fliplr(p), origin='lower')
        # plt.show()
        gamma_options = {
            'lower_percent_dose_cutoff': int(cuoff),
            'interp_fraction': int(interfra),  # Should be 10 or more for more accurate results
            'random_subset': None,
            'max_gamma': float(maxgamma),
            'local_gamma': False,
            'quiet': True
        }
        if moreinfo :
            gamma_options['quiet']=False
        local_gamma = False
        max_dose = 0
        if 'local' in maxdose:
            gamma_options['local_gamma'] = True
        elif 'glo' in maxdose:
            pass
        else:
            gamma_options['global_normalisation '] = float(maxdose)

        criterias = dosediscrit.split(',')
        gammalist=[]
        criterialist=[]
        for criteria in criterias:
            dosecrit = criteria.split('/')[0]
            discrit = criteria.split('/')[1]
            criterialist.append(dosecrit+'%/'+discrit+'mm')
            gamma_options['dose_percent_threshold'] = int(dosecrit)
            gamma_options['distance_mm_threshold'] = int(discrit)

            #start_time = time.time()
            gamma = pymedphys.gamma(
                axes_reference, dose_reference,
                axes_evaluation, dose_evaluation_fx,
                **gamma_options)

            valid_gamma = gamma[~np.isnan(gamma)]
            gammavalue=len(valid_gamma[valid_gamma <= 1]) / len(valid_gamma) * 100
            print(
                f"Criteria {dosecrit}%/{discrit}mm Passing Rate(\u03B3<=1): {len(valid_gamma[valid_gamma <= 1]) / len(valid_gamma) * 100}%")
            #print("cputime ", time.time() - start_time)
            gammalist.append(round(gammavalue,2))
        # write files
        No_firstline = True
        try:
            with open(saveresultas, 'r') as file_read:
                if 'ref' in file_read.readlines()[0]:
                    No_firstline = False
        except:
            pass
        with open(saveresultas, 'a+') as file_save:
            # file_save.writelines(str(datetime.today()) + ' ' + str(datetime.utcnow()) + '\n')
            if (No_firstline):
                file_save.writelines('addtionalinfo   reference   compare  criteria Passing-rate\n')
            additionalinfo_replace=''
            if additionalinfo!=None:
                additionalinfo_replace=additionalinfo.replace('%',' ')
            file_save.writelines(additionalinfo_replace+' '+dose1[-9:]+' '+dose2[-9:]+' ')
            for temp in range(0,len(gammalist)):
                file_save.writelines(str(criterialist[temp])+' '+str(gammalist[temp])+ '% ')
            file_save.write("\n")

    def fun_1Dgamma(self):
        reference = np.genfromtxt('dose_film_1D_z0mm.csv', delimiter=',', skip_header=1)
        evaluation = np.genfromtxt('dose_MC_1D_z0mm.csv', delimiter=',', skip_header=1)

        axis_reference = reference[:, 0]  # 1st column is x in mm
        dose_reference = reference[:, 1]  # 2nd column is dose in Gy/MU

        axis_evaluation = evaluation[:, 0]
        dose_evaluation = evaluation[:, 1]
        gamma_options = {
            'dose_percent_threshold': 1,
            'distance_mm_threshold': 1,
            'lower_percent_dose_cutoff': 10,
            'interp_fraction': 10,  # Should be 10 or more for more accurate results
            'max_gamma': 2,
            'random_subset': None,
            'local_gamma': False,  # False indicates global gamma is calculated
            'ram_available': 2 ** 29  # 1/2 GB
        }
        # for global dose normalization, the maximum reference dose is used
        # but in TG218, it said usually the prescribed dose or the maximum dose in a plan (evaluation) is used
        gamma = pymedphys.gamma(
            axis_reference, dose_reference,
            axis_evaluation, dose_evaluation,
            **gamma_options)
        print(gamma)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--ref", required=True,
                        help="reference dose file (nrrd/DCM)")
    parser.add_argument("-c", "--comp", required=True,
                        help="compare dose file (nrrd/DCM)")
    parser.add_argument("-dodi", "--dosediscrit", required=False, nargs='?',
                        help="criteria for dose deviation and distance to agreement, e.g. 3/3,2/2,1/1. default:3/3", default=3/3)
    parser.add_argument("-cu", "--cutoff", required=False, nargs='?',
                        help="The percent lower dose cutoff above which gamma will be calculated. Only applied to the reference grid. Default: 10",
                        default=10)
    parser.add_argument("-m", "--maxdose", required=False, nargs='?',
                        help="max dose for percentage calculation, one of 'local','global',or number,default:global", default='global')
    parser.add_argument("-i", "--interfra", required=False, nargs='?',
                        help=" The fraction which gamma distance threshold is divided into for interpolation. Defaults to 5 (0.6mm step for 3mm criteria), less=faster",
                        default=5)
    parser.add_argument("-mg", "--maxgamma", required=False, nargs='?',
                        help="Once a search distance is reached that would give gamma values larger than this parameter, the search stops. less=faster. >1, default=1.1",
                        default=1.1)
    parser.add_argument("-f", "--fraction", required=False, nargs='?',
                        help="No. of fractions for compared dose, scale by multiply, default:1", default=1)
    parser.add_argument("-s", "--saveas", required=False, nargs='?',
                        help="Save the gamma result to file the path and name of the file. Default=./gammaresults.txt", default='./gammaresults.txt')
    parser.add_argument("-p", "--prone", required=False, action='store_true',
                        help="prone case, reference nrrd will be flipped lr and un,default:not active",default=False)
    parser.add_argument("-a", "--addition", required=False, nargs='?',
                        help="set additional information for output such as beam name, plan name, patient ID etc.")
    parser.add_argument("-info", "--info", required=False, action='store_true',
                        help="active more information mode", default=False)
    args = parser.parse_args()
    print('start a new analysis')
    gammaana=class_gammaanalysis()
    gammaana.fun_gamma_analysis(args.addition,args.ref, args.comp, args.dosediscrit, args.cutoff, args.maxdose, args.interfra, args.maxgamma,
                       args.fraction, args.saveas, args.prone, args.info)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
