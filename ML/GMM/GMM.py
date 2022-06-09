from scipy.stats import multivariate_normal
import numpy as np 

class Gaussian: 
    def __init__(self, mean = np.ones((2)), covariance = np.eye(2), dim = 2):
        self.mean, self.covariance = mean, covariance

    def __call__(self, x):
        return multivariate_normal.pdf(x, mean=self.mean, cov=self.covariance)
    
class GMM():
    def __init__(self, k, dim):
        init_mu = [np.random.random((dim))-0.5 for i in range(k)]
        init_sigma = [np.eye(dim) for i in range(k)]
        self.Gaussians = [Gaussian(init_mu[i], init_sigma[i], dim) for i in range(k)]
        self.pi = np.ones((k)) / k
        self.data = []
        self.K = k 
        self.N = 0
        self.dim = dim
    
    def add_data(self, x):
        self.data.append(x)
        self.N += 1
    
    def prob(self, x_idx, g_idx):
        return self.Gaussians[g_idx](self.data[x_idx])

    def E_step(self):
        z_new = np.zeros((self.N, self.K))
        for i, data in enumerate(self.data):
            for j, g in enumerate(self.Gaussians):
                z_new[i, j] = g(data) * self.pi[j]
            z_new[i] /= z_new[i].sum()
        self.z = z_new
        print(self.z)

    def M_step(self):
        self.pi = self.z.sum(axis=0)
        self.pi = self.pi / self.pi.sum()
        for j in range(self.K):
            mean_vector = np.zeros((self.dim))
            for i in range(self.N):
                mean_vector += self.z[i, j] * self.data[i]
            mean_vector /= self.z[:, j].sum()
            self.Gaussians[j].mean = mean_vector
            cov_mat = np.zeros((self.dim, self.dim))
            for i in range(self.N):
                cov_mat += self.z[i, j] * np.outer((self.data[i] - mean_vector), (self.data[i] - mean_vector))
            cov_mat /= self.z[:, j].sum()
            self.Gaussians[j].covariance = cov_mat
    
    def log_likelihood(self):
        return round(sum([np.log(sum([g(x) * self.pi[j] for (j, g) in enumerate(self.Gaussians)])) for (i, x) in enumerate(self.data)]), 2)
        
    def print_result(self):
        for g in self.Gaussians:
            print(f"Gaussian mu = {g.mean.round(2)}\n         S = {np.array2string(g.covariance.round(2), prefix='         S = ')}")
        print()


if __name__ == '__main__':
    gmm = GMM(2, 1)
    gmm.Gaussians[0].mean = np.float32([1])
    gmm.Gaussians[0].covariance = np.float32([1])
    gmm.Gaussians[1].mean = np.float32([10])
    gmm.Gaussians[1].covariance = np.float32([1])

    X = [1, 2, 3, 7, 9, 10]
    for i in range(len(X)):
        gmm.add_data(X[i])

    print(f"Likelihood = {gmm.log_likelihood()}")
    gmm.print_result()
    for steps in range(2):
        gmm.E_step()
        gmm.M_step()
        print(f"Likelihood = {gmm.log_likelihood()}")
        gmm.print_result()
        