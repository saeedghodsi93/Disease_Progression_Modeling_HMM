
# import the modules
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.special import rel_entr, kl_div
import pickle


# store the raw data in file
def dthmm_store_raw(I, J, L, N_vals, M_vals, H, z, y, u):
  
  # put the results in a dict and store them in a file
  dict_raw = {"I": I, "J": J, "L": L, "N_vals": N_vals, "M_vals": M_vals, "H": H, "z": z, "y": y, "u": u}
  with open('DTHMM/Data/raw.pickle', 'wb') as handle:
    pickle.dump(dict_raw, handle, protocol=pickle.HIGHEST_PROTOCOL)
  
  return


# load the raw data from file
def dthmm_load_raw():
  
  # load data from file
  with open('DTHMM/Data/raw.pickle', 'rb') as handle:
    dict_raw = pickle.load(handle)
  
  return dict_raw


# store the EM results in file
def dthmm_store_results(pi, Q, mu, eta, pi_hat_vals, Q_hat_vals, mu_hat_vals, eta_hat_vals):
  
  # put the results in a dict and store them in a file
  dict_res = {"pi": pi, "Q": Q, "mu": mu, "eta": eta, "pi_hat_vals": pi_hat_vals, "Q_hat_vals": Q_hat_vals, "mu_hat_vals": mu_hat_vals, "eta_hat_vals": eta_hat_vals}
  with open('DTHMM/Data/res.pickle', 'wb') as handle:
    pickle.dump(dict_res, handle, protocol=pickle.HIGHEST_PROTOCOL)
  
  return


# load the stored EM results from file
def dthmm_load_results():
  
  # load data from file
  with open('DTHMM/Data/res.pickle', 'rb') as handle:
    dict_res = pickle.load(handle)
  
  return dict_res


# get the data from dthmm.R and visualize it
def dthmm_visualize_data(H, z, y, u):
  
  # choose a patient for visualization
  patient_idx = 7
  
  # extract the relevant data
  patient_H = round(H[patient_idx])
  patient_z = z[patient_idx,:]
  patient_y = y[patient_idx,:]
  patient_u = u[patient_idx,:]
  patient_z = [int(element) for element in patient_z if ~np.isnan(element)]
  patient_y = [int(element) for element in patient_y if ~np.isnan(element)]
  patient_u = [int(element) for element in patient_u if ~np.isnan(element)]
  
  fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(10,10), gridspec_kw={'height_ratios': [1, 2, 1]})
  l1 = ax1.step(range(1,patient_H+1), patient_z, linestyle='-.', marker='o', markersize=3, color="maroon", label="z", where="post")
  l2 = ax2.step(range(1,patient_H+1), patient_y, linestyle='-', marker='o', markersize=3, color="peru", label="y", where="post")
  l3 = ax3.step(range(1,patient_H+1), patient_u, linestyle='--', marker='o', markersize=3, color="seagreen", label="u", where="post")
  ax1.get_xaxis().set_ticks([])
  ax2.get_xaxis().set_ticks([])
  ax3.get_xaxis().set_ticks(range(1,patient_H+1))
  ax1.get_yaxis().set_ticks([1,2,3])
  ax2.get_yaxis().set_ticks([0,2,4,6,8])
  ax3.get_yaxis().set_ticks([0,1,2])
  ax1.set_ylim(1-0.2, 3.2)
  ax2.set_ylim(-0.5, 9.5)
  ax3.set_ylim(-0.2, 2.2)
  ax1.set_xlim(0.8, patient_H+0.2)
  ax2.set_xlim(0.8, patient_H+0.2)
  ax3.set_xlim(0.8, patient_H+0.2)
  for t in range(1,patient_H+1):
    ax1.axvline(x=t, color='gray', linestyle='--', alpha=0.5)
    ax2.axvline(x=t, color='gray', linestyle='--', alpha=0.5)
    ax3.axvline(x=t, color='gray', linestyle='--', alpha=0.5)
  ax1.grid(axis='y')
  ax2.grid(axis='y')
  ax3.grid(axis='y')
  fig.suptitle("The underlying health state (z), physician observation (y), and intervention (u) variables", fontsize = "medium")
  ax1.set(ylabel="z")
  ax2.set(ylabel="y")
  ax3.set(ylabel="u")
  fig.supxlabel("time", fontsize = "medium")
  plt.subplots_adjust(top=0.93)
  plt.show()
  
  fig.savefig("DTHMM/Results/synthetic_data.pdf", bbox_inches='tight')
  plt.close(fig)
  
  return


# KL divergence
def dthmm_kl_divergence(p, q):
  
  dist = sum(kl_div(p, q))
  
  return dist


# calculate the distances between the true variables and the em predicted variables
def dthmm_calc_distances(I, J, L, N_vals, M_vals, num_iterations, pi, Q, mu, eta, pi_hat_vals, Q_hat_vals, mu_hat_vals, eta_hat_vals):
  
  # calculate the distances for each iteration
  dist_N = []
  dist_M = []
  dist_itr = []
  dist_pi = []
  dist_Q = []
  dist_mu = []
  dist_eta = []
  for N_itr in range(len(N_vals)):
    for M_itr in range(len(M_vals)):
      for em_itr in range(num_iterations):
        
        # add the value of N, M, and em iteration 
        dist_N.append(N_vals[N_itr])
        dist_M.append(M_vals[M_itr])
        dist_itr.append(em_itr)
      
        # extract the parameters corresponding to each iteration
        pi_hat_itr = pi_hat_vals[N_itr, M_itr, em_itr,]
        Q_hat_itr = Q_hat_vals[N_itr, M_itr, em_itr,]
        mu_hat_itr = mu_hat_vals[N_itr, M_itr, em_itr]
        eta_hat_itr = eta_hat_vals[N_itr, M_itr, em_itr,]
        
        # calculate the distance for pi
        dist_pi.append(np.sqrt(np.mean(np.square(pi - pi_hat_itr))))
      
        # calculate the distance for Q
        temp = 0
        for l in range(L):
          for i in range(I):
            for k in range(I):
              temp += np.float_power(Q[i,k,l] - Q_hat_itr[i,k,l], 2)
        dist_Q.append(np.sqrt(temp / (I * I * L)))
        
        # calculate the distance for mu
        dist_mu.append(np.sqrt(np.mean(np.square(mu - mu_hat_itr))))
        
        # calculate the distance for eta
        dist_eta.append(np.sqrt(np.mean(np.square(eta - eta_hat_itr))))
    
  # create and return a data frame
  dict_dist = {"N": dist_N, "M": dist_M, "itr": dist_itr, "pi": dist_pi, "Q": dist_Q, "mu": dist_mu, "eta": dist_eta}
  df_dist = pd.DataFrame(dict_dist)

  return df_dist


# get the results from dthmm.R and visualize them
def dthmm_visualize_results(I, J, L, N_vals, M_vals, pi, Q, mu, eta, pi_hat_vals, Q_hat_vals, mu_hat_vals, eta_hat_vals):
  
  # convert the lists to numpy array and correct the integers
  N_vals = [round(N) for N in N_vals]
  M_vals = [round(M) for M in M_vals]
  pi_hat_vals = np.asarray((pi_hat_vals))
  Q_hat_vals = np.asarray((Q_hat_vals))
  mu_hat_vals = np.asarray((mu_hat_vals))
  eta_hat_vals = np.asarray((eta_hat_vals))
  
  # the number of N values and iterations
  num_iterations = pi_hat_vals.shape[2]
  
  # calculate the distances
  df_dist = dthmm_calc_distances(round(I), round(J), round(L), N_vals, M_vals, num_iterations, pi, Q, mu, eta, pi_hat_vals, Q_hat_vals, mu_hat_vals, eta_hat_vals)
    
  # define the lists of colors, linestyles, and legends
  color_vals = [["tan", "maroon"], ["purple", "cadetblue"]]
  linestyle_vals = [["-", "-"], ["--", "--"]]
  marker_vals = [["o", "o"], ["x", "x"]]
  legend_vals = []
  for N in N_vals:
    for M in M_vals:
      if M == 0:
        legend_vals.append("N="+str(N)+", direct")
      else:
        legend_vals.append("N="+str(N)+", M="+str(M))
        
  # plot the distances for pi
  fig = plt.figure(figsize=(5,4))
  for N_itr in range(len(N_vals)):
    for M_itr in range(len(M_vals)):
      df_dist_N_M = df_dist[(df_dist["N"]==N_vals[N_itr]) & (df_dist["M"]==M_vals[M_itr])]
      plt.plot(np.asarray(df_dist_N_M["itr"]), np.asarray(df_dist_N_M["pi"]), linestyle=linestyle_vals[N_itr][M_itr], marker=marker_vals[N_itr][M_itr], markersize=3, color=color_vals[N_itr][M_itr])
  plt.xlim(left=-1, right=num_iterations)
  plt.ylim(bottom=-0.005)
  plt.grid(True, linestyle='--', linewidth=1, alpha=0.7)
  plt.legend(labels=legend_vals, loc="upper right")
  plt.title("Convergence of pi to the true value")
  plt.xlabel("iteration")
  plt.ylabel("RMSE (pi)")
  plt.show()
  fig.savefig("DTHMM/Results/convergence_pi.pdf", bbox_inches='tight')
  plt.close(fig)
  
  fig = plt.figure(figsize=(5,4))
  for N_itr in range(len(N_vals)):
    for M_itr in range(len(M_vals)):
      df_dist_N_M = df_dist[(df_dist["N"]==N_vals[N_itr]) & (df_dist["M"]==M_vals[M_itr])]
      plt.plot(np.asarray(df_dist_N_M["itr"]), np.asarray(df_dist_N_M["Q"]), linestyle=linestyle_vals[N_itr][M_itr], marker=marker_vals[N_itr][M_itr], markersize=3, color=color_vals[N_itr][M_itr])
  plt.xlim(left=-1, right=num_iterations)
  plt.ylim(bottom=-0.002)
  plt.grid(True, linestyle='--', linewidth=1, alpha=0.7)
  plt.legend(labels=legend_vals, loc="upper right")
  plt.title("Convergence of Q to the true value")
  plt.xlabel("iteration")
  plt.ylabel("RMSE (Q)")
  plt.show()
  fig.savefig("DTHMM/Results/convergence_Q.pdf", bbox_inches='tight')
  plt.close(fig)
  
  fig = plt.figure(figsize=(5,4))
  for N_itr in range(len(N_vals)):
    for M_itr in range(len(M_vals)):
      df_dist_N_M = df_dist[(df_dist["N"]==N_vals[N_itr]) & (df_dist["M"]==M_vals[M_itr])]
      plt.plot(np.asarray(df_dist_N_M["itr"]), np.asarray(df_dist_N_M["mu"]), linestyle=linestyle_vals[N_itr][M_itr], marker=marker_vals[N_itr][M_itr], markersize=3, color=color_vals[N_itr][M_itr])
  plt.xlim(left=-1, right=num_iterations)
  plt.ylim(bottom=-0.005)
  plt.grid(True, linestyle='--', linewidth=1, alpha=0.7)
  plt.legend(labels=legend_vals, loc="upper right")
  plt.title("Convergence of mu to the true value")
  plt.xlabel("iteration")
  plt.ylabel("RMSE (mu)")
  plt.show()
  fig.savefig("DTHMM/Results/convergence_mu.pdf", bbox_inches='tight')
  plt.close(fig)
  
  fig = plt.figure(figsize=(5,4))
  for N_itr in range(len(N_vals)):
    for M_itr in range(len(M_vals)):
      df_dist_N_M = df_dist[(df_dist["N"]==N_vals[N_itr]) & (df_dist["M"]==M_vals[M_itr])]
      plt.plot(np.asarray(df_dist_N_M["itr"]), np.asarray(df_dist_N_M["eta"]), linestyle=linestyle_vals[N_itr][M_itr], marker=marker_vals[N_itr][M_itr], markersize=3, color=color_vals[N_itr][M_itr])
  plt.xlim(left=-1, right=num_iterations)
  plt.ylim(bottom=-0.003)
  plt.grid(True, linestyle='--', linewidth=1, alpha=0.7)
  plt.legend(labels=legend_vals, loc="upper right")
  plt.title("Convergence of eta to the true value")
  plt.xlabel("iteration")
  plt.ylabel("RMSE (eta)")
  plt.show()
  fig.savefig("DTHMM/Results/convergence_eta.pdf", bbox_inches='tight')
  plt.close(fig)
  
  return

