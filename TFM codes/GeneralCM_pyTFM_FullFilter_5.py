# THIS code is to be used for Prostate Cancer Readouts


from pyTFM.TFM_functions import calculate_deformation, TFM_tractions, strain_energy_points, contractillity
from pyTFM.plotting import show_quiver, plot_continuous_boundary_stresses
from pyTFM.stress_functions import lineTension
from pyTFM.grid_setup_solids_py import interpolation, prepare_forces, grid_setup, FEM_simulation, find_borders
from pyTFM.utilities_TFM import round_flexible
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.morphology import binary_fill_holes
import os
import pandas as pd
import xlsxwriter 
import scipy.io as sio
import math

#import sys
#sys.path.append("D:\RAW DATA\EC_TFM_3kPa_PC3wwoAR_R1881_LNCaPwwoR1881_CSS_CTOrange_48h_220519")


# Import filter for noise removal 

#from low_passfilter import low_passfilter

# Import Low Pass Filter

## calculating a deformation field
## calculating a traction forces

# important parameters:
ps1 =  0.5687  # 0.5687  0.284 pixel size of tehe image of the beads Previous * 0.284  # NucImage = 0.0901876
im1_shape = (472,472) # dimensions of the image of the beads # Image will be 512x512 
young = 3000 # Young's modulus of the substrate
sigma = 0.49 # Poisson's ratio of the substrate
h = 100 # height of the substrate in µm, "infinite" is also accepted
window_size_var = 50#80   
overlap_var =  47#6
step_size = window_size_var - overlap_var
std_factor_var = 5

# Initialize Results Folder
results = {"strain_energy": [],  "rmst": [],"average_displacement": [], "contractile_force": [], "ps1": [], "ps2":[], "im1_shape": [] , "h": [], "young":[], "window_size":[], "overlap":[], "step_size":[], "pixel_to_array_dims":[], "PIV_output_size":[], "pixel_start_PIV_buffer":[],  "trim_buffer_win":[], "final_pixel_start":[] }

results["ps1"] = ps1
results["im1_shape"]  = im1_shape
results["h"] = h
results["young"] = young
results["window_size"] = window_size_var
results["overlap"] = overlap_var
results["step_size_ppw"] = step_size 
ps2 = ps1*step_size;
results["ps2"] = ps2
results["std_factor"] =std_factor_var
pixel_to_array_dims = int(im1_shape[0]/step_size)
results["pixel_to_array_dims"] = pixel_to_array_dims


#full_frame_background; 
if  np.mod(window_size_var, step_size):
    pixel_start_PIV_buffer = math.floor(window_size_var/2) 
else :
    pixel_start_PIV_buffer = math.floor(window_size_var/2)+1
        
results["pixel_start_PIV_buffer"] = pixel_start_PIV_buffer 
    

  #"ps2":[], , "h": [], "young":[], "window_size":[], "overlap":[], "trim_factor":[], "window_buffer":[]}

##############################################################################################################################################################################################################

#create a loop to go over files from time-lapse and append strain energy values
# pip install numpy==1.19.5 --user
# G:\My Drive\Master's\Research\Prostate\pyTFM\BB_LNCaP_CSS_FBS_ML_3Day_inc_FN_TC_3kPa_DEC12_2020\LNCaP 3kPa CSS
os.chdir(r"D:\RAW DATA\EH_OMTCTFM_PC3AR_3kPa_10nMR1881_mOrange_48h_HEAT_06102022\PC3AR_10nMR1881_3kPa")
folder = os.getcwd()
print(folder)

#initialize the dataframes for every position in between 
         
 # Initiate the lists 
rmst_list = []
strain_energy_list = []
average_displacement_list = []
deformationU_df_collection_list = []
deformationV_df_collection_list = []
tractionU_df_collection_list = []
tractionV_df_collection_list = []
max_disp_check_list = []

#Loop through positions of interest 

use_filter = 1; 
series_index = [1,2,3,4,5,6,7,8,9,10]
for series_range in range(len(series_index)): #6
 
    im_path1 = os.path.join(folder, "series" + "_" + str(series_index[series_range]) + "/" + "Series"+str(series_index[series_range]) ,"Null.tif")
    #print(im_path1)
    
    a= len(os.listdir(os.path.join(folder, "series" + "_" + str(series_index[series_range]) + "/" + "Series"+str(series_index[series_range])+"/Beads/")))
    rmst_df = pd.DataFrame(columns = ['Timepoint','Position'+str(series_index[series_range])]) 
    strain_energy_df = pd.DataFrame(columns = ['Timepoint','Position'+str(series_index[series_range])])
    average_displacement_df = pd.DataFrame(columns = ['Timepoint','Position'+str(series_index[series_range])])
    contractile_force_df = pd.DataFrame(columns = ['Timepoint','Position'+str(series_index[series_range])])
    max_disp_check_df = pd.DataFrame(columns = ['Timepoint','Position'+str(series_index[series_range])])
    
    tractionU_df_collection = {}
    tractionV_df_collection = {}
    
    frame_lim = 11
    mylist = [1]
    for x in range(1, a+1): #mylist: range(1, a+1) # range(1, a+1):#frame_lim+1): # range(1, a+1) mylist : 
        
        im_path2 = os.path.join(folder,"series" + "_" + str(series_index[series_range]) + "/" + "Series"+str(series_index[series_range])+"/Beads/","Beads_Forced_"+str(x).zfill(4)+ ".tif")
        print(im_path2)
     
    #use window_size = 100, overlap = 80 for 1024x1024 images or window_size = 50, overlap = 40 for 1024x1024 images         
        
        u,v, mask_val, mask_std = calculate_deformation(im_path1, im_path2,  window_size = window_size_var, overlap = overlap_var, std_factor = std_factor_var)      
        PIV_buffer = pixel_to_array_dims-np.array(u.shape); 
        PIV_output_size = u.shape[0]# OLD ps2 = ps1 * np.mean(np.array(im1_shape) / np.array(u.shape)) # pixel size of of the deformation field
        results["PIV_output_size"] = PIV_output_size
        
        #anchor_point_pixel = window_size/step-size
       
                # Create Displacement folder  
        if not os.path.exists(os.path.join(folder,"series" + "_" + str(series_index[series_range]) + "/" +"Series"+str(series_index[series_range])+"/Displacement")):
            os.mkdir(os.path.join(folder,"series" + "_" + str(series_index[series_range]) + "/" +"Series"+str(series_index[series_range])+"/Displacement"))
        # Create Uncorrected Displacement Folder 
        if not os.path.exists(os.path.join(folder,"series" + "_" + str(series_index[series_range]) + "/" +"Series"+str(series_index[series_range])+"/Displacement/d_arrays_uncorrected")):
            os.mkdir(os.path.join(folder,"series" + "_" + str(series_index[series_range]) + "/" +"Series"+str(series_index[series_range])+"/Displacement/d_arrays_uncorrected"))       
        #Create Corrected Displalcement FOlder  
        if not os.path.exists(os.path.join(folder,"series" + "_" + str(series_index[series_range]) + "/" +"Series"+str(series_index[series_range])+"/Displacement/d_arrays_corrected")):
            os.mkdir(os.path.join(folder,"series" + "_" + str(series_index[series_range]) + "/" +"Series"+str(series_index[series_range])+"/Displacement/d_arrays_corrected"))

#########################      
#FILTER        # Create Create FFT Folder 
        if not os.path.exists(os.path.join(folder,"series" + "_" + str(series_index[series_range]) + "/" +"Series"+str(series_index[series_range])+"/Displacement/FFT")):
            os.mkdir(os.path.join(folder,"series" + "_" + str(series_index[series_range]) + "/" +"Series"+str(series_index[series_range])+"/Displacement/FFT"))
        
        filepath = os.path.join(folder,"series" + "_" + str(series_index[series_range]) + "/" + "Series"+str(series_index[series_range])+"/Displacement/FFT") 
        filename ="Displacement"+str(x).zfill(3)
       
       # if in use, uncomment and comment accordingly 
       # u_filter, v_filter = low_passfilter(u, v, filepath, filename, x)
  ##############################      
        
       # Factor is the pixel size
       # Save the non-cropped/filtered displacement map
        fig1, ax = show_quiver(u*ps1, v*ps1, cbar_str=" ", alpha=1, plot_cbar=True, scale_ratio=0.1 ,filter=[0, 5], width=0.003)# plotting
        fig1.savefig(os.path.join(folder,"series" + "_" + str(series_index[series_range]) + "/" +"Series"+str(series_index[series_range])+"/Displacement/d_arrays_uncorrected/raw_Displacement"+str(x).zfill(3)+".tif"))
        
     
        t_dim = u.shape[0]
        t_reduce = math.floor(int(t_dim*0.1))
        u_trim = u[(t_reduce-1):(t_dim-1-t_reduce), (t_reduce-1):(t_dim-1-t_reduce)]
        v_trim = v[(t_reduce-1):(t_dim-1-t_reduce), (t_reduce-1):(t_dim-1-t_reduce)]
        
       # trim_buffer = np.shape(u) - np.shape(u_trim) 
        
        # Plot and save the trimmed matrix 
        fig1b, ax = show_quiver(u_trim*ps1, v_trim*ps1,  cbar_str=" ", vmin = 0, vmax= 10, alpha=1, plot_cbar=True, scale_ratio=0.1, filter=[0, 5], width=0.003)# plotting
        fig1b.savefig(os.path.join(folder,"series" + "_" + str(series_index[series_range]) + "/" +"Series"+str(series_index[series_range])+"/Displacement/trim_Displacement"+str(x).zfill(3)+".tif"))
        
        #Save trimmed matrices 
        displacement_mat =  {'u':u*ps1, 'v':v*ps1}
        displacement_mat_crop = {'u_filt':u_trim*ps1, 'v_filt':v_trim*ps1}
        sio.savemat(os.path.join(folder, "series" + "_" + str(series_index[series_range]) + "/" +"Series"+str(series_index[series_range])+ "/Displacement/d_arrays_uncorrected/raw_displacement_mat" + str(x).zfill(3) + ".mat"), displacement_mat)
        sio.savemat(os.path.join(folder, "series" + "_" + str(series_index[series_range]) + "/" +"Series"+str(series_index[series_range])+ "/Displacement/d_arrays_corrected/trim_displacement_mat" + str(x).zfill(3) + ".mat"), displacement_mat_crop)
        
    #loading a pseudo-mask that defines the area used for measuring the force generation
        mask = np.ones(im1_shape)
        mask = interpolation(mask, dims=u.shape)

        # Calculate Tractions from the original displacement matrix 
        tx, ty = TFM_tractions(u, v, pixelsize1=ps1, pixelsize2=ps2, h=h, young=young, sigma=sigma)
        contractile_force, proj_x, proj_y, center = contractillity(tx, ty, ps2, mask) # 2.01*10**-6 N 
        
       
        t_dim = tx.shape[0]
        t_reduce = math.floor(int(t_dim*0.1))
        tx_trim = tx[(t_reduce-1):(t_dim-1-t_reduce), (t_reduce-1):(t_dim-1-t_reduce)]
        ty_trim = ty[(t_reduce-1):(t_dim-1-t_reduce), (t_reduce-1):(t_dim-1-t_reduce)]
        
        trim_buffer = t_reduce # = u.shape[0] - u_trim.shape[0];
      
        results["trim_buffer_win"] = trim_buffer
        results["final_pixel_start"] = pixel_start_PIV_buffer+trim_buffer*step_size
        
       # For now, report only the traction
        traction_mat = {'tzx':tx, 'tzy':ty}
        traction_mat_trim = {'tzx_filt':tx_trim, 'tzy':ty_trim}        

        ##########Filtered 
        #tx_filt, ty_filt = TFM_tractions(u_filter, v_filter, pixelsize1=ps1, pixelsize2=ps2, h=h, young=young, sigma=sigma)
        #contractile_force_filt, proj_x_filt, proj_y_filt, center_filt = contractillity(tx, ty, ps2, mask) # 2.01*10**-6 N 
        
        #t_dim = tx_filt.shape[0]
        #t_reduce = round(int(t_dim*0.1))
        #tx_trim_filt = tx_filt[(t_reduce-1):(t_dim-1-t_reduce), (t_reduce-1):(t_dim-1-t_reduce)]
        #ty_trim_filt = ty_filt[(t_reduce-1):(t_dim-1-t_reduce), (t_reduce-1):(t_dim-1-t_reduce)]

#########################################
       
        fig2, ax = show_quiver(tx, ty, cbar_str="", alpha=0.85,plot_cbar=True, scale_ratio=0.1,filter=[0, 5], width=0.003)
        
        fig2b, ax = show_quiver(tx_trim, ty_trim, cbar_str="", vmin=0, vmax=500, alpha=0.85, plot_cbar=True, scale_ratio=0.1,filter=[0, 5], width=0.003)
       # t_dim = tx.shape[0]
       # t_reduce = round(int(t_dim*0.1, decimals = none)
        
        tractionU_df_collection[x] = pd.DataFrame(data=tx[1:,1:], index=tx[1:,0], columns=tx[0,1:])
        tractionV_df_collection[x] = pd.DataFrame(data=ty[1:,1:], index=ty[1:,0], columns=ty[0,1:])
        

        if not os.path.exists(os.path.join(folder,"series" + "_" + str(series_index[series_range]) + "/" +"Series"+str(series_index[series_range])+"/Traction")):
            os.mkdir(os.path.join(folder,"series" + "_" + str(series_index[series_range]) + "/" +"Series"+str(series_index[series_range])+"/Traction"))
            
        if not os.path.exists(os.path.join(folder,"series" + "_" + str(series_index[series_range]) + "/" +"Series"+str(series_index[series_range])+"/Traction/t_arrays_uncorrected")):
            os.mkdir(os.path.join(folder,"series" + "_" + str(series_index[series_range]) + "/" +"Series"+str(series_index[series_range])+"/Traction/t_arrays_uncorrected"))
            
        if not os.path.exists(os.path.join(folder,"series" + "_" + str(series_index[series_range]) + "/" +"Series"+str(series_index[series_range])+"/Traction/t_arrays_corrected")):
            os.mkdir(os.path.join(folder,"series" + "_" + str(series_index[series_range]) + "/" +"Series"+str(series_index[series_range])+"/Traction/t_arrays_corrected"))       
        
        fig2.savefig(os.path.join(folder,"series" + "_" + str(series_index[series_range]) + "/" +"Series"+str(series_index[series_range])+"/Traction/t_arrays_uncorrected/raw_traction"+str(x).zfill(3)+".tif"))        
        sio.savemat(os.path.join(folder, "series" + "_" + str(series_index[series_range]) + "/" +"Series"+str(series_index[series_range])+ "/Traction/t_arrays_uncorrected/raw_traction_mat" + str(x).zfill(3) + ".mat"), traction_mat)
        
        fig2b.savefig(os.path.join(folder,"series" + "_" + str(series_index[series_range]) + "/" +"Series"+str(series_index[series_range])+"/Traction/trim_traction"+str(x).zfill(3)+".tif"))        
        sio.savemat(os.path.join(folder,"series" + "_" + str(series_index[series_range]) + "/" + "Series"+str(series_index[series_range])+ "/Traction/t_arrays_corrected/trim_traction_mat" + str(x).zfill(3) + ".mat"), traction_mat_trim)
        
    
        #mask = plt.imread(os.path.join(folder,"mask-colony","monolayer-"+str(x+1)+".png")).astype(bool)
    	#mask = binary_fill_holes(mask) # the mask should be a single patch without holes
    	#mask = interpolation(mask, dims=u_trim.shape)

        
        ########### COMPUTE SUMMARY STATISTICS 
        u = u_trim
        v = v_trim 
        tx = tx_trim 
        ty = ty_trim 
        
        energy_points = strain_energy_points(u, v, tx, ty, ps1, ps2) # J/pixel
    	#strain_energy = np.sum(energy_points[mask]) # 1.92*10**-13 J
        strain_energy = np.sum(energy_points) # 1.92*10**-13 J --> Clayton mod
        rmst = np.sqrt(np.mean(tx**2+ty**2))
        average_displacement = np.sqrt(np.mean(u**2+v**2))*ps1
        max_disp_check = np.max(np.sqrt(u**2+v**2)*ps1)
        #results["strain_energy"].append(strain_energy)
        #results["rmst"].append(rmst)
        #results["contractile_force"].append(contractile_force)
        #results["average_displacement"].append(average_displacement)
        
        
        strain_energy_df.loc[len(strain_energy_df)] = [x, strain_energy]
        contractile_force_df.loc[len(contractile_force_df)] = [x, contractile_force]
        rmst_df.loc[len(rmst_df)] = [x, rmst]
        average_displacement_df.loc[len(average_displacement_df)] =  [x, average_displacement]
        max_disp_check_df.loc[len(max_disp_check_df)] = [x, max_disp_check]
    
    rmst_list.append(rmst_df.drop('Timepoint', axis=1))
    strain_energy_list.append(strain_energy_df.drop('Timepoint', axis=1))
    average_displacement_list.append(average_displacement_df.drop('Timepoint', axis=1))  
    max_disp_check_list.append(max_disp_check_df.drop('Timepoint', axis=1))
    
rmst_result_df = pd.concat(rmst_list, axis=1)
rmst_result_array = rmst_result_df.to_numpy()
#rmst_result_array = np.array(rmst_result_df.to_numpy(), dtype=object)
strain_energy_result_df = pd.concat(strain_energy_list, axis=1)
strain_energy_result_array =strain_energy_result_df.to_numpy()

max_disp_check_result_df =  pd.concat(max_disp_check_list, axis=1)
max_disp_check_result_array =  max_disp_check_result_df.to_numpy()

average_displacement_result_df = pd.concat(average_displacement_list, axis =1)
average_displacement_result_array = average_displacement_result_df.to_numpy()

results["rmst"] = rmst_result_array #test #rmst_result_df
results["strain_energy"] =strain_energy_result_array
results["average_displacement"] = average_displacement_result_array
results["max_disp_check"] = max_disp_check_result_array

sio.savemat(os.path.join(folder, "results.mat"), results)


print("Hello")