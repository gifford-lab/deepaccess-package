import os
import keras
import pickle
import numpy as np
from keras.models import load_model, Sequential, Model
from keras.layers import Dense, Input, Add
from keras import optimizers
from deepaccess.ensemble_utils import *
from deepaccess.train.CNN import *
from deepaccess.interpret.importance_utils import saliency

class DeepAccessModel:
    def __init__(self, outdir):
        self.model_paths = []
        self.accuracies = {}
        self.outdir = outdir
        ensure_dir(outdir)
        self.ensemble = None

    def _make_ensemble(self):
        models = [load_model(mp + "/model.h5") for mp in self.model_paths]
        model_input = Input(shape=models[0].input.shape[1:])

        new_models = []
        for model in models:
            l = model_input
            for i in range(0, len(model.layers)):
                l = model.layers[i](l)
            new_models.append(Model(model_input, l))

        acc = [self.accuracies[mp.split("/")[-1]] for mp in self.model_paths]
        acc_norm = np.array(acc) / sum(acc)
        outputs = [model.outputs[0] * acc_norm[mi] for mi, model in enumerate(new_models)]
        y = Add()(outputs)
        model = Model(model_input, y, name="DeepAccessEnsemble")
        return model

    def load(self, model_folder):
        self.model_paths = sorted(
            [
                model_folder + "/" + d
                for d in os.listdir(model_folder)
                if os.path.isdir(model_folder + "/" + d)
            ],
            key=lambda x: int(x.split("_")[-1]),
        )
        with open(model_folder + "/model_acc.pkl", "rb") as f:
            initial_accuracies = pickle.load(f)
        for key in list(initial_accuracies.keys()):
            self.accuracies[key.split("/")[-1]] = initial_accuracies[key]
        self.nclasses = np.loadtxt(model_folder + "/test.txt").shape[1]
        self.ensemble = self._make_ensemble()

    def predict(self, X):
        return self.ensemble.predict(X)

    def saliency_input(self, X, c1, c2, n=5, batch_size=256):
        saliencyX = np.zeros((X.shape[0], X.shape[1]))
        for mi, model in enumerate(self.model_paths):
            saliencyX += np.sum(
                saliency(
                    0,
                    model + "/model.h5",
                    X,
                    c1,
                    c2,
                    n=n,
                    lambda_reg=0.1,
                    batch_size=batch_size,
                )
                * X,
                axis=2,
            )
            saliencyX *= self.accuracies[model.split("/")[-1]]
        return saliencyX

    def train(self, X, y, n_epochs=10, verbose=0):
        models = [['conv','globalpool'],
                  ['conv','conv','globalpool'],
                  ['conv','conv','globalpool'],
                  ['conv','conv','globalpool'],
                  ['conv','conv','globalpool'],
                  ['conv','conv','globalpool'],
                  ['conv','conv','globalpool'],
                  ['conv','conv','globalpool'],
                  ['conv','conv','globalpool'],
                  ['conv','conv','globalpool']
        ]
        filters =[None,
                  (200,8),
                  (200,16),
                  (200,32),
                  (400,8),
                  (400,16),
                  (400,32),
                  (600,8),
                  (600,16),
                  (600,32)]
        
        sample_weights = np.ones((X.shape[0],))
        accuracies = {}
        self.nclasses = y.shape[1]
        trained_models = []
        for mi, model in enumerate(models):
            if filters[mi] != None:
                new_cnn = CNN(model,
                              X.shape[1:],
                              y.shape[1],
                              conv_filter_number=filters[mi][0],
                              conv_filter_size=filters[mi][1])
            else:
                new_cnn = CNN(model,
                              X.shape[1:],
                              y.shape[1])
            
            history = new_cnn.train(X, y, sample_weights, n_epochs=n_epochs)
            model_folder = self.outdir + "/" + "_".join(model)+"_"+str(mi)
            ensure_dir(model_folder)
        
            loss, train_acc = new_cnn.evaluate(
                X, y, batch_size=250, verbose=verbose
            )
            
            accuracies[model_folder] = train_acc
            new_cnn.save(model_folder + "/model.h5")
            trained_models.append(model_folder)
            np.save(model_folder + "/history.npy", history.history)
            np.save(model_folder + "/sample_weights.npy", sample_weights)
            with open(model_folder + "/model_summary.txt", "w") as f:
                f.write(new_cnn.to_yaml())
                f.write("\nTraining Accuracy: " + str(train_acc) + "\n")

            sample_weights += np.linalg.norm(
                y - new_cnn.predict(X, batch_size=250)
            )
            del new_cnn
        self.model_paths = trained_models
        for key in list(accuracies.keys()):
            self.accuracies[key.split("/")[-1]] = accuracies[key]
        self.ensemble = self._make_ensemble()

        with open(self.outdir + "/model_acc.pkl", "wb") as f:
            pickle.dump(accuracies, f)
