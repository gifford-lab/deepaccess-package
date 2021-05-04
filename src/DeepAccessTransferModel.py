import os
import keras
import pickle
from ensemble_utils import ensure_dir
import numpy as np
from keras.models import load_model, Sequential, Model
from keras.layers import Dense, Input, Add
from keras import optimizers
from importance_utils import saliency


class DeepAccessTransferModel:
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
        outputs = [
            model.outputs[0] * acc_norm[mi]
            for mi, model in enumerate(new_models)
        ]
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

    def retrain(self, X, y, n_epochs=1, verbose=0):
        sample_weights = np.ones((X.shape[0],))
        accuracies = {}
        self.nclasses = y.shape[1]
        retrained_models = []
        for mi, model in enumerate(self.model_paths):
            cnn = keras.models.load_model(model + "/model.h5")
            # add up to last layer
            new_cnn = Sequential()
            for layer in cnn.layers[:-1]:
                new_cnn.add(layer)
                # add new layer with new output
            new_cnn.add(Dense(y.shape[1], activation="sigmoid"))
            adam = optimizers.Adam(
                lr=1e-4,
                clipnorm=0.5,
                decay=(1e-4 / 100.0))
            new_cnn.compile(
                optimizer=adam,
                loss="binary_crossentropy",
                metrics=["accuracy"]
            )

            model_folder = self.outdir + "/" + model.split("/")[-1]
            ensure_dir(model_folder)
            callbacks = [
                keras.callbacks.EarlyStopping(monitor="val_loss", patience=3),
                keras.callbacks.History(),
            ]
            history = new_cnn.fit(
                x=X,
                y=y,
                epochs=n_epochs,
                shuffle=True,
                validation_split=0.2,
                batch_size=250,
                verbose=verbose,
                callbacks=callbacks,
            )
            loss, train_acc = new_cnn.evaluate(
                X,
                y,
                batch_size=250,
                verbose=verbose)
            accuracies[model_folder] = train_acc
            new_cnn.save(model_folder + "/model.h5")
            retrained_models.append(model_folder)
            np.save(model_folder + "/history.npy", history.history)
            np.save(model_folder + "/sample_weights.npy", sample_weights)
            with open(model_folder + "/model_summary.txt", "w") as f:
                f.write(new_cnn.to_yaml())
                f.write("\nTraining Accuracy: " + str(train_acc) + "\n")

            sample_weights += np.linalg.norm(
                y - new_cnn.predict(X, batch_size=250)
            )

            del cnn
            del new_cnn
        self.model_paths = retrained_models
        for key in list(accuracies.keys()):
            self.accuracies[key.split("/")[-1]] = accuracies[key]
        self.ensemble = self._make_ensemble()

        with open(self.outdir + "/model_acc.pkl", "wb") as f:
            pickle.dump(accuracies, f)
