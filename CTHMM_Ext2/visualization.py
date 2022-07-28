
# import the modules
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.special import rel_entr, kl_div
import pickle


# store the raw data in file
def cthmm_ext2_store_raw(I, J, L, N_vals, H, a, O, tau_true, tau_obs, z_true, z_obs, z_acc, y_obs, u_obs):
  
  # put the results in a dict and store them in a file
  dict_raw = {"I": I, "J": J, "L": L, "N_vals": N_vals, "H": H, "a": a, "O": O, "tau_true": tau_true, "tau_obs": tau_obs, "z_true": z_true, "z_obs": z_obs, "z_acc": z_acc, "y_obs": y_obs, "u_obs": u_obs}
  with open('cthmm_ext2/Data/raw.pickle', 'wb') as handle:
    pickle.dump(dict_raw, handle, protocol=pickle.HIGHEST_PROTOCOL)
  
  return


# load the raw data from file
def cthmm_ext2_load_raw():
  
  # load data from file
  with open('cthmm_ext2/Data/raw.pickle', 'rb') as handle:
    dict_raw = pickle.load(handle)
  
  return dict_raw


# store the EM results in file
def cthmm_ext2_store_results(pi, Q, mu, eta, eta_prime, pi_hat_vals, Q_hat_vals, mu_hat_vals, eta_hat_vals, eta_prime_hat_vals):
  
  # put the results in a dict and store them in a file
  dict_res = {"pi": pi, "Q": Q, "mu": mu, "eta": eta, "eta_prime": eta_prime, "pi_hat_vals": pi_hat_vals, "Q_hat_vals": Q_hat_vals, "mu_hat_vals": mu_hat_vals, "eta_hat_vals": eta_hat_vals, "eta_prime_hat_vals": eta_prime_hat_vals}
  with open('cthmm_ext2/Data/res.pickle', 'wb') as handle:
    pickle.dump(dict_res, handle, protocol=pickle.HIGHEST_PROTOCOL)
  
  return


# load the stored EM results from file
def cthmm_ext2_load_results():
  
  # load data from file
  with open('cthmm_ext2/Data/res.pickle', 'rb') as handle:
    dict_res = pickle.load(handle)
  
  return dict_res


# get the data from cthmm_ext2.R and visualize it
def cthmm_ext2_visualize_data(H, a, O, tau_true, tau_obs, z_true, z_obs, z_acc, y_obs, u_obs):
  
  # choose a patient for visualization
  patient_idx = 7
  
  # extract the relevant data
  patient_tau_true = tau_true[patient_idx,:]
  patient_z_true = z_true[patient_idx,:]
  patient_tau_obs = tau_obs[patient_idx,:]
  patient_z_obs = z_obs[patient_idx,:]
  patient_z_acc = z_acc[patient_idx,:]
  patient_y_obs = y_obs[patient_idx,:]
  patient_u_obs = u_obs[patient_idx,:]
  patient_tau_true = [element if ~np.isnan(element) else element for element in patient_tau_true]
  patient_z_true = [int(element) if ~np.isnan(element) else element for element in patient_z_true]
  patient_tau_obs = [element if ~np.isnan(element) else element for element in patient_tau_obs]
  patient_z_obs = [int(element) if ~np.isnan(element) else element for element in patient_z_obs]
  patient_z_acc = [int(element) if ~np.isnan(element) else element for element in patient_z_acc]
  patient_y_obs = [int(element) if ~np.isnan(element) else element for element in patient_y_obs]
  patient_u_obs = [int(element) if ~np.isnan(element) else element for element in patient_u_obs]
  
  fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, figsize=(10,12), gridspec_kw={'height_ratios': [1, 1, 2, 1]})
  l1 = ax1.step(patient_tau_true, patient_z_true, linestyle='--', marker='o', markersize=3, color="maroon", label="z", where="post")
  l2 = ax2.step(patient_tau_obs, patient_z_acc, linestyle='--', marker='o', markersize=3, color="darkmagenta", label="z.acc", where="post")
  l3 = ax3.step(patient_tau_obs, patient_y_obs, linestyle='--', marker='o', markersize=3, color="peru", label="y", where="post")
  l4 = ax4.step(patient_tau_obs, patient_u_obs, linestyle='--', marker='o', markersize=3, color="seagreen", label="u", where="post")
  ax1.get_xaxis().set_ticks([])
  ax2.get_xaxis().set_ticks([])
  ax3.get_xaxis().set_ticks([])
  ax1.get_yaxis().set_ticks([1,2,3])
  ax2.get_yaxis().set_ticks([1,2,3])
  ax3.get_yaxis().set_ticks([0,2,4,6,8])
  ax4.get_yaxis().set_ticks([0,1,2])
  ax1.set_ylim(1-0.2, 3.2)
  ax2.set_ylim(1-0.2, 3.2)
  ax3.set_ylim(-0.5, 9.5)
  ax4.set_ylim(-0.2, 2.2)
  ax1.set_xlim(-2, np.nanmax(patient_tau_obs)+2)
  ax2.set_xlim(-2, np.nanmax(patient_tau_obs)+2)
  ax3.set_xlim(-2, np.nanmax(patient_tau_obs)+2)
  ax4.set_xlim(-2, np.nanmax(patient_tau_obs)+2)
  for t in patient_tau_obs:
    ax1.axvline(x=t, color='gray', linestyle='--', alpha=0.5)
    ax2.axvline(x=t, color='gray', linestyle='--', alpha=0.5)
    ax3.axvline(x=t, color='gray', linestyle='--', alpha=0.5)
    ax4.axvline(x=t, color='gray', linestyle='--', alpha=0.5)
  ax1.grid(axis='y')
  ax2.grid(axis='y')
  ax3.grid(axis='y')
  ax4.grid(axis='y')
  fig.suptitle("The underlying health state (z), physician observation (y), and intervention (u) variables", fontsize = "medium")
  ax1.set(ylabel="z")
  ax2.set(ylabel="z.acc")
  ax3.set(ylabel="y")
  ax4.set(ylabel="u")
  fig.supxlabel("time", fontsize = "medium")
  plt.subplots_adjust(top=0.93)
  plt.show()
  
  fig.savefig("cthmm_ext2/Results/synthetic_data.pdf", bbox_inches='tight')
  plt.clf()
  
  return


# KL divergence
def cthmm_ext2_kl_divergence(p, q):
  
  dist = sum(kl_div(p, q))
  
  return dist


# calculate the distances between the true variables and the em predicted variables
def cthmm_ext2_calc_distances(I, J, L, N_vals, num_iterations, pi, Q, mu, eta, eta_prime, pi_hat_vals, Q_hat_vals, mu_hat_vals, eta_hat_vals, eta_prime_hat_vals):
  
  # calculate the distances for each iteration
  dist_N = []
  dist_itr = []
  dist_pi = []
  dist_Q = []
  dist_mu = []
  dist_eta = []
  dist_eta_prime = []
  for N_itr in range(len(N_vals)):
    for em_itr in range(num_iterations):
      
      # add the value of N and em iteration 
      dist_N.append(N_vals[N_itr])
      dist_itr.append(em_itr)
    
      # extract the parameters corresponding to each iteration
      pi_hat_itr = pi_hat_vals[N_itr, em_itr,]
      Q_hat_itr = Q_hat_vals[N_itr, em_itr,]
      mu_hat_itr = mu_hat_vals[N_itr, em_itr]
      eta_hat_itr = eta_hat_vals[N_itr, em_itr,]
      eta_prime_hat_itr = eta_prime_hat_vals[N_itr, em_itr,]
      
      # calculate the distance for pi
      dist_pi.append(cthmm_ext2_kl_divergence(pi, pi_hat_itr))
    
      # calculate the distance for Q
      temp = 0
      for l in range(L):
        for i in range(I):
          for k in range(I):
            if k != i:
              temp += np.float_power(Q[i,k,l] - Q_hat_itr[i,k,l], 2)
      dist_Q.append(np.sqrt(temp / (I * (I-1) * L)))
      
      # calculate the distance for mu
      dist_mu.append(np.sqrt(np.mean(np.square(mu - mu_hat_itr))))
      
      # calculate the distance for eta
      dist_eta.append(np.sqrt(np.mean(np.square(eta - eta_hat_itr))))
      
      # calculate the distance for eta_prime
      dist_eta_prime.append(np.sqrt(np.mean(np.square(eta_prime - eta_prime_hat_itr))))
    
  # create and return a data frame
  dict_dist = {"N": dist_N, "itr": dist_itr, "pi": dist_pi, "Q": dist_Q, "mu": dist_mu, "eta": dist_eta, "eta_prime": dist_eta_prime}
  df_dist = pd.DataFrame(dict_dist)

  return df_dist


# get the results from cthmm_ext2.R and visualize them
def cthmm_ext2_visualize_results(I, J, L, N_vals, pi, Q, mu, eta, eta_prime, pi_hat_vals, Q_hat_vals, mu_hat_vals, eta_hat_vals, eta_prime_hat_vals):
  
  # convert the lists to numpy array and correct the integers
  N_vals = [round(N) for N in N_vals]
  pi_hat_vals = np.asarray((pi_hat_vals))
  Q_hat_vals = np.asarray((Q_hat_vals))
  mu_hat_vals = np.asarray((mu_hat_vals))
  eta_hat_vals = np.asarray((eta_hat_vals))
  eta_prime_hat_vals = np.asarray((eta_prime_hat_vals))
  
  # the number of N values and iterations
  num_iterations = pi_hat_vals.shape[1]
  
  # calculate the distances
  df_dist = cthmm_ext2_calc_distances(round(I), round(J), round(L), N_vals, num_iterations, pi, Q, mu, eta, eta_prime, pi_hat_vals, Q_hat_vals, mu_hat_vals, eta_hat_vals, eta_prime_hat_vals)
  
  # define the lists of colors, linestyles, and legends
  color_vals = ["tan", "maroon", "seagreen", "cadetblue"]
  linestyle_vals = ["-", "-.", "--", ":"]
  legend_vals = ["N="+str(N) for N in N_vals]
  
  # plot the distances for pi
  fig = plt.figure(figsize=(5,4))
  for N_itr in range(len(N_vals)):
    df_dist_N = df_dist[df_dist["N"]==N_vals[N_itr]]
    plt.plot(np.asarray(df_dist_N["itr"]), np.asarray(df_dist_N["pi"]), linestyle=linestyle_vals[N_itr], marker='o', markersize=3, color=color_vals[N_itr])
  plt.xlim(left=-1, right=num_iterations)
  plt.ylim(bottom=-0.005)
  plt.grid(True, linestyle='--', linewidth=1, alpha=0.7)
  plt.legend(labels=legend_vals, loc="upper right")
  plt.title("Convergence of pi to the true value")
  plt.xlabel("iteration")
  plt.ylabel("RMSE (pi)")
  plt.show()
  fig.savefig("cthmm_ext2/Results/convergence_pi.pdf", bbox_inches='tight')
  plt.clf()
  
  fig = plt.figure(figsize=(5,4))
  for N_itr in range(len(N_vals)):
    df_dist_N = df_dist[df_dist["N"]==N_vals[N_itr]]
    plt.plot(np.asarray(df_dist_N["itr"]), np.asarray(df_dist_N["Q"]), linestyle=linestyle_vals[N_itr], marker='o', markersize=3, color=color_vals[N_itr])
  plt.xlim(left=-1, right=num_iterations)
  plt.ylim(bottom=-0.002)
  plt.grid(True, linestyle='--', linewidth=1, alpha=0.7)
  plt.legend(labels=legend_vals, loc="upper left")
  plt.title("Convergence of delta to the true value")
  plt.xlabel("iteration")
  plt.ylabel("RMSE (delta)")
  plt.show()
  fig.savefig("cthmm_ext2/Results/convergence_delta.pdf", bbox_inches='tight')
  plt.clf()
  
  fig = plt.figure(figsize=(5,4))
  for N_itr in range(len(N_vals)):
    df_dist_N = df_dist[df_dist["N"]==N_vals[N_itr]]
    plt.plot(np.asarray(df_dist_N["itr"]), np.asarray(df_dist_N["mu"]), linestyle=linestyle_vals[N_itr], marker='o', markersize=3, color=color_vals[N_itr])
  plt.xlim(left=-1, right=num_iterations)
  plt.ylim(bottom=-0.005)
  plt.grid(True, linestyle='--', linewidth=1, alpha=0.7)
  plt.legend(labels=legend_vals, loc="upper right")
  plt.title("Convergence of mu to the true value")
  plt.xlabel("iteration")
  plt.ylabel("RMSE (mu)")
  plt.show()
  fig.savefig("cthmm_ext2/Results/convergence_mu.pdf", bbox_inches='tight')
  plt.clf()
  
  fig = plt.figure(figsize=(5,4))
  for N_itr in range(len(N_vals)):
    df_dist_N = df_dist[df_dist["N"]==N_vals[N_itr]]
    plt.plot(np.asarray(df_dist_N["itr"]), np.asarray(df_dist_N["eta"]), linestyle=linestyle_vals[N_itr], marker='o', markersize=3, color=color_vals[N_itr])
  plt.xlim(left=-1, right=num_iterations)
  plt.ylim(bottom=-0.003)
  plt.grid(True, linestyle='--', linewidth=1, alpha=0.7)
  plt.legend(labels=legend_vals, loc="upper right")
  plt.title("Convergence of eta to the true value")
  plt.xlabel("iteration")
  plt.ylabel("RMSE (eta)")
  plt.show()
  fig.savefig("cthmm_ext2/Results/convergence_eta.pdf", bbox_inches='tight')
  plt.clf()
  
  fig = plt.figure(figsize=(5,4))
  for N_itr in range(len(N_vals)):
    df_dist_N = df_dist[df_dist["N"]==N_vals[N_itr]]
    plt.plot(np.asarray(df_dist_N["itr"]), np.asarray(df_dist_N["eta_prime"]), linestyle=linestyle_vals[N_itr], marker='o', markersize=3, color=color_vals[N_itr])
  plt.xlim(left=-1, right=num_iterations)
  plt.ylim(bottom=-0.003)
  plt.grid(True, linestyle='--', linewidth=1, alpha=0.7)
  plt.legend(labels=legend_vals, loc="upper right")
  plt.title("Convergence of eta.prime to the true value")
  plt.xlabel("iteration")
  plt.ylabel("RMSE (eta.prime)")
  plt.show()
  fig.savefig("cthmm_ext2/Results/convergence_eta_prime.pdf", bbox_inches='tight')
  plt.clf()
  
  return

