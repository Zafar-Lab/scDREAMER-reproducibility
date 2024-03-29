import tensorflow as tf2

import tensorflow.compat.v1 as tf
tf.disable_v2_behavior() 

import scanpy as sc
import pandas as pd
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import LabelEncoder
import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics.cluster import normalized_mutual_info_score as nmi

def read_h5ad(data_path, batch, cell_type, name, sparseIP, hvg=2000):
    print('updated hvg')
    Ann = sc.read_h5ad(data_path)
    Ann.layers["counts"] = Ann.X.copy()

    sc.pp.normalize_total(Ann, target_sum=1e4)
    sc.pp.log1p(Ann)
    Ann.raw = Ann 
    
    sc.pp.highly_variable_genes(
        Ann, 
        flavor="seurat", 
        n_top_genes=hvg,
        batch_key=batch,
        subset=True)
  
    if (sparseIP == 1):
        df_final = pd.DataFrame.sparse.from_spmatrix(Ann.X) # Lung, Simulation1, Simulation 2
    else:
        df_final = pd.DataFrame(Ann.X) # Immune , Pan, Healthy Heart

    df_final = df_final.reset_index(drop = True)
    data = df_final.to_numpy()
    
    if cell_type:
        labels = Ann.obs[cell_type].to_list()

    #AJ: Convert to categorical instead of this...
    t_ = Ann.obs[batch] #.to_list()
    batch_info = np.array([[i] for i in t_]) # for other datasets

    #batch_info = np.array(Ann.obs[batch].astype("category").reset_index(drop = True)).reshape(-1,1) 
    enc = OneHotEncoder(handle_unknown='ignore')
    #batch_info_enc = enc.fit_transform(batch_info).toarray()
  
    enc.fit(batch_info.reshape(-1, 1))
    batch_info_enc = enc.transform(batch_info.reshape(-1, 1)).toarray()
    
    return data, labels, batch_info_enc, batch_info
"""
def read_h5ad(data_path, B, C):
    Ann = sc.read_h5ad(data_path)
    Ann.layers["counts"] = Ann.X.copy()

    sc.pp.normalize_total(Ann, target_sum=1e4)
    sc.pp.log1p(Ann)
    Ann.raw = Ann 
    
    sc.pp.highly_variable_genes(
        Ann, 
        flavor="seurat", 
        n_top_genes=hvg,
        batch_key=B,
        subset=True)
  
    df_final = pd.DataFrame.sparse.from_spmatrix(Ann.X) # Lung, Simulation1, Simulation 2
    #df_final = pd.DataFrame(Ann.X) # Immune , Pan

    df_final = df_final.reset_index(drop = True)
    data = df_final.to_numpy()
    labels = Ann.obs[C].to_list()

    #AJ: Convert to categorical instead of this...
    t_ = Ann.obs[B] #.to_list()
    batch_info = np.array([[i] for i in t_]) # for other datasets

    #batch_info = np.array(Ann.obs[B].astype("category").reset_index(drop = True)).reshape(-1,1) 
    enc = OneHotEncoder(handle_unknown='ignore')
    #batch_info_enc = enc.fit_transform(batch_info).toarray()
  
    enc.fit(batch_info.reshape(-1, 1))
    batch_info_enc = enc.transform(batch_info.reshape(-1, 1)).toarray()

    return data, labels, batch_info_enc, batch_info
"""
# Leaky Relu
def lrelu(x, alpha = 0.2, name='lrelu'):
    return tf.maximum(x, alpha*x)

def dense(x, inp_dim, out_dim, name = 'dense'):

    with tf.variable_scope(name, reuse=None): # earlier only tf
        weights = tf.get_variable("weights", shape=[inp_dim, out_dim],
                                  initializer = #tf2.keras.initializers.glorot_uniform(seed = 0))
                                  tf2.initializers.GlorotUniform()) # contrib: tf.contrib.layers.xavier_initializer()
        
        bias = tf.get_variable("bias", shape=[out_dim], initializer = tf.constant_initializer(0.0))
        
        # initializer= tf2.initializers.GlorotUniform(); same as Xavier's initializer; tf.contrib.layers.xavier_initializer()    
        out = tf.add(tf.matmul(x, weights), bias, name='matmul')
        return out

def load_gene_mtx(dataset_name, name, transform = True, count = True, actv = 'sig', batch = "batch", \
                  cell_type = "cell_type", sparseIP = 0):
    print('came in load_gene')
    #B = "tech"
    #C = "celltype"
    data, labels, batch_info_enc, batch_info = read_h5ad(dataset_name, batch, cell_type, name, sparseIP)
         
    if count == False:
        data = np.log2(data + 1) #np.ones(data.shape[0], data.shape[1]))
        #data = np.log2(data)
        
        if actv == 'lin':
            scale = 1.0
        else:
            scale = np.max(data)
        data = data / scale           

    ord_enc = LabelEncoder()
    labels  = ord_enc.fit_transform(labels)
    print ('here', labels)

    unique, counts = np.unique(labels, return_counts = True)
    dict(zip(unique, counts))
    
    total_size = data.shape[0]

    if count == False:
        return data, data, scale, labels, labels, batch_info_enc, batch_info_enc, batch_info

    return data, data, labels, labels, labels

"""
def load_gene_mtx(dataset_name, transform = True, count = True, actv = 'sig'):

    data, labels, batch_info_enc, batch_info = read_h5ad(path + data_name, B, C)
        
    print ('Shape of data is: ', data.shape)

    if transform:
        data = transform_01_gene_input(data)
        print('Data Transformed, entries in [0, 1] !'.format(data_path))
    else:
        if count == False:
            data = np.log2(data+1)

            if actv == 'lin':
                scale = 1.0
            else:
                scale = np.max(data)
            data = data / scale           

    ord_enc = LabelEncoder()
    labels  = ord_enc.fit_transform(labels)
    print ('here', labels)

    unique, counts = np.unique(labels, return_counts = True)
    dict(zip(unique, counts))
    
    total_size = data.shape[0]

    if count == False:
        return data, data, scale, labels, labels, batch_info_enc, batch_info_enc, batch_info

    return data, data, labels, labels, labels
    
"""

# AJ: 20 may

# Zero-inflated negative binomial (ZINB) model is for modeling count variables with excessive zeros and it is usually for overdispersed count outcome variables.

def zinb_model(self, x, mean, inverse_dispersion, logit, eps=1e-4): 

    # 1e8 should be of same dimensions as other parameters....                 
    expr_non_zero = - tf.nn.softplus(- logit) \
                    + tf.log(inverse_dispersion + eps) * inverse_dispersion \
                    - tf.log(inverse_dispersion + mean + eps) * inverse_dispersion \
                    - x * tf.log(inverse_dispersion + mean + eps) \
                    + x * tf.log(mean + eps) \
                    - tf.lgamma(x + 1) \
                    + tf.lgamma(x + inverse_dispersion) \
                    - tf.lgamma(inverse_dispersion) \
                    - logit 
    
    expr_zero = - tf.nn.softplus( - logit) \
                + tf.nn.softplus(- logit + tf.log(inverse_dispersion + eps) * inverse_dispersion \
                                  - tf.log(inverse_dispersion + mean + eps) * inverse_dispersion) 

    template = tf.cast(tf.less(x, eps), tf.float32)
    expr =  tf.multiply(template, expr_zero) + tf.multiply(1 - template, expr_non_zero)
    return tf.reduce_sum(expr, axis=-1)

def eval_cluster_on_test_(self, epoch):

    # Embedding points in the test data to the latent space
    inp_encoder = self.data_test
    labels = self.labels_test
    batch_label = self.batch_test
            
    #latent_matrix = self.sess.run(self.z, feed_dict = {self.x_input: inp_encoder, self.batch_input: batch_label, self.keep_prob: 1.0})
    start = 0
    end = 15000 # size of each pass
    latent_matrix = self.sess.run(self.z, feed_dict = {self.x_input: inp_encoder[start:end], self.batch_input: batch_label[start:end], self.keep_prob: 1.0})

    while (end < len(inp_encoder)):

      #print ("hi")
      start = end
      end = min(end + 15000, len(inp_encoder))

      mat = self.sess.run(self.z, feed_dict = {self.x_input: inp_encoder[start:end], self.batch_input: batch_label[start:end], self.keep_prob: 1.0})
      latent_matrix = np.concatenate((latent_matrix, mat), axis = 0)
    
    print ('latent_matrix shape', latent_matrix.shape)
    print (labels.shape)
    
    Ann = sc.AnnData(inp_encoder)
    Ann.obsm['final_embeddings'] = latent_matrix
    Ann.obs['group'] = labels.astype(str)
    
    sc.pp.neighbors(Ann, use_rep = 'final_embeddings') #use_rep = 'final_embeddings'
    sc.tl.umap(Ann)
    img = sc.pl.umap(Ann, color = 'group', frameon = False) # cells
    print(img)
    
    np.savetxt(self.name + 'latent_matrix_' + str(epoch) +'.csv', latent_matrix, delimiter=",")
    
    Ann.obs['batch'] = self.batch_info.astype(str)
    img2 = sc.pl.umap(Ann, color = 'batch', frameon = False)
    print(img2)

    K = np.size(np.unique(labels))   
    kmeans = KMeans(n_clusters=K, random_state=0).fit(latent_matrix)
    y_pred = kmeans.labels_

    print('Computing NMI ...')
    NMI = nmi(labels.flatten(), y_pred.flatten())
    print('Done !')

    print('NMI = {}'. 
          format(NMI)) 

def eval_cluster_on_test(self, epoch):

    # Embedding points in the test data to the latent space
    inp_encoder = self.data_test
    labels = self.labels_test
    batch_label = self.batch_test
            
    #latent_matrix = self.sess.run(self.z, feed_dict = {self.x_input: inp_encoder, self.batch_input: batch_label, self.keep_prob: 1.0})
    
    start = 0
    end = 15000 # size of each pass
    latent_matrix = self.sess.run(self.z, feed_dict = {self.x_input: inp_encoder[start:end], self.batch_input: batch_label[start:end], self.keep_prob: 1.0})

    while (end < len(inp_encoder)):

      #print ("hi")
      start = end
      end = min(end + 15000, len(inp_encoder))

      mat = self.sess.run(self.z, feed_dict = {self.x_input: inp_encoder[start:end], self.batch_input: batch_label[start:end], self.keep_prob: 1.0})
      latent_matrix = np.concatenate((latent_matrix, mat), axis = 0)
    
    print ('latent_matrix shape', latent_matrix.shape)
    print (labels.shape)
    
    Ann = sc.AnnData(inp_encoder)
    Ann.obsm['final_embeddings'] = latent_matrix
    Ann.obs['group'] = labels.astype(str)
    
    sc.pp.neighbors(Ann, use_rep = 'final_embeddings') #use_rep = 'final_embeddings'
    sc.tl.umap(Ann)
    img = sc.pl.umap(Ann, color = 'group', frameon = False) # cells
    print(img)
    
    np.savetxt(self.name + 'latent_matrix_' + str(epoch) +'.csv', latent_matrix, delimiter=",")
    
    Ann.obs['batch'] = self.batch_info.astype(str)
    img2 = sc.pl.umap(Ann, color = 'batch', frameon = False)
    print(img2)

    K = np.size(np.unique(labels))   
    kmeans = KMeans(n_clusters=K, random_state=0).fit(latent_matrix)
    y_pred = kmeans.labels_

    print('Computing NMI ...')
    NMI = nmi(labels.flatten(), y_pred.flatten())
    print('Done !')

    print('NMI = {}'. 
          format(NMI)) 

