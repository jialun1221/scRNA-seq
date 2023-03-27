import scanpy
from sklearn.model_selection import train_test_split
import seaborn as sns
from sklearn.metrics import ConfusionMatrixDisplay
from sklearn import metrics
from sklearn import preprocessing
import sklearn.metrics as metrics
import pandas as pd

from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report

import torch
import numpy as np
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import torchvision
from torch.utils.data import TensorDataset

use_gpu = torch.cuda.is_available()
device = torch.device("cuda:0" if use_gpu else "cpu")

from tqdm import tqdm
from numpy import mean
from numpy import std
from sklearn.datasets import make_classification
from sklearn.model_selection import cross_val_score

import matplotlib
from matplotlib import pyplot as plt

def logistic_regression(adata, seed):
    if 'level_0' in adata.obs.columns:
        print('no need to reset')
    else:
        adata.obs = adata.obs.reset_index()

    if 'X_pca' in adata.obsm:
        X = adata.obsm['X_pca'].X
    else:
        X = adata.X
    #print(adata.obsm['X_pca'])
    #X = adata.obsm['X_pca'].X
    #print(X)
    #print(type(X))
    #X = np.array(X)
    y = adata.obs['disease__ontology_label'].replace({"normal": "0", "Parkinson disease": "1"})
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=0)

    use_gpu = torch.cuda.is_available()
    device = torch.device("cuda:0" if use_gpu else "cpu")
    np.random.seed(seed)  # Set the random seed of numpy for the data split.
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)

    model = LogisticRegression()
    # print(X_train.shape,y_train.shape)
    # print(type(X_train))
    # print(type(y_train))
    model.fit(X_train, y_train)
    prediction_test = model.predict(X_test)

    accuracy = metrics.accuracy_score(y_test, prediction_test)
    print("Logistic Regression Accuracy = ", metrics.accuracy_score(y_test, prediction_test))
    print("End Logistic Regression")

    #plot
    confusion_matrix = metrics.confusion_matrix(y_test, prediction_test)
    cmd = ConfusionMatrixDisplay(confusion_matrix, display_labels=['Healthy', 'PD'])
    cmd.plot(cmap="crest", colorbar=False)

    fig = cmd.ax_.get_figure()
    fig.set_figwidth(3)
    fig.set_figheight(3)

    plt.rcParams["figure.dpi"] = 500
    plt.rcParams["axes.labelcolor"] = 'black'  # now useless

    plt.rc('font', size=16)
    # Set the axes title font size
    plt.rc('axes', titlesize=10)
    # Set the axes labels font size
    plt.rc('axes', labelsize=10)
    # Set the font size for x tick labels
    plt.rc('xtick', labelsize=10)
    # Set the font size for y tick labels
    plt.rc('ytick', labelsize=10)

    plt.show()
    # return prediction_test, accuracy
    return accuracy

#prediction_test, accuracy = logistic_regression()

# def plot_logistic_regression(y_test, prediction_test, adata):
    # X = pd.DataFrame(adata.obsm['X_pca'])
    # y = adata.obs['disease__ontology_label'].replace({"normal": "0", "Parkinson disease": "1"})
    # X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=0)

    # confusion_matrix = metrics.confusion_matrix(y_test, prediction_test)
    # cmd = ConfusionMatrixDisplay(confusion_matrix, display_labels=['Healthy', 'PD'])
    # cmd.plot(cmap="crest")
    #
    # fig = cmd.ax_.get_figure()
    # fig.set_figwidth(3)
    # fig.set_figheight(3)
    #
    # plt.rcParams["figure.dpi"] = 500
    # plt.rcParams["axes.labelcolor"] = 'black'  # now useless

def random_forest(adata, seed):
    print('start of random forest')
    if 'level_0' in adata.obs.columns:
        print('no need to reset')
    else:
        adata.obs = adata.obs.reset_index()

    if 'X_pca' in adata.obsm:
        X = adata.obsm['X_pca'].X
    else:
        X = adata.X
    y = adata.obs['disease__ontology_label'].replace({"normal": "0", "Parkinson disease": "1"})

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=0)

    np.random.seed(seed)  # Set the random seed of numpy for the data split.
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)

    forest = RandomForestClassifier(random_state=0, n_estimators=100)
    forest.fit(X_train, y_train)

    prediction_test = forest.predict(X_test)
    accuracy = accuracy_score(y_test, prediction_test)
    print("Accuracy is:" + str(accuracy))

    #plot
    print("reaching plot")
    confusion_matrix = metrics.confusion_matrix(y_test, prediction_test)
    cmd = ConfusionMatrixDisplay(confusion_matrix, display_labels=['Healthy', 'PD'])
    cmd.plot(cmap="Blues", colorbar=False)

    fig = cmd.ax_.get_figure()
    fig.set_figwidth(3)
    fig.set_figheight(3)

    # fig.set_ylabel(fontsize = 10)
    plt.rcParams["figure.dpi"] = 500
    plt.rcParams["axes.labelcolor"] = 'black'  # now useless
    plt.show()
    plt.savefig('RF_confusion_matrix.pdf')
    # return prediction_test, accuracy
    return accuracy


def DeepNeuralNet(adata, seed):
    print('Deep Neural Net')
    if 'level_0' in adata.obs.columns:
        print('no need to reset')
    else:
        adata.obs = adata.obs.reset_index()

    if 'X_pca' in adata.obsm:
        #print('reached here')
        X = adata.obsm['X_pca'].X
    else:
        print('reached here istead!')
        X = adata.X
    #print('X_train ', X.shape)
    y = adata.obs['disease__ontology_label'].replace({"normal": "0", "Parkinson disease": "1"})

    np.random.seed(seed)  # Set the random seed of numpy for the data split.
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)

    # train test split
    from sklearn.model_selection import train_test_split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=0)
    #set different seeds for randomized split of training and testing (validation)
    #test should be the same
    #cross validation randomly does it for you
    #validation should be smallest
    #plot variations (boxplot), and dot for each accuracy of 5 runs
    print('X_train size: ', X_train.shape, 'X: ', X_train.shape)

    #print(type(X_train))
    #print(type(X_test))
    X_train = pd.DataFrame(X_train)
    X_test = pd.DataFrame(X_test)
    y_train = pd.DataFrame(y_train)
    y_test = pd.DataFrame(y_test)
    # print(type(X_train))
    # print(type(X_test))

    # Take some data from test set for validation
    # Get the validation dataset splitted
    X_val = X_train.head(992)
    X_train = X_train.tail(5000)
    y_val = y_train.head(992)
    y_train = y_train.tail(5000)
    print('shape summary: ', X_val.shape, X_test.shape, y_val.shape, y_test.shape, X_train.shape, y_train.shape)
    print(type(X_val), type(X_test), type(y_val), type(y_test), type(X_train), type(y_train)) #THey are all dataframe

    # Type conversion
    X_val = X_val.astype(np.float64)
    y_val = y_val.astype(np.float64)
    X_val = X_val.to_numpy()
    y_val = y_val.to_numpy()
    # print("here" + str(y_val.shape))

    # Note: SoftMax requires a 2x1 matrix for the label. We turn the above `y_val`, `y_train` and `y_test` numpy array into a matrix.
    #y_val = np.stack([y_val, 1 - y_val])
    #y_val = y_val.transpose()
    y_val = y_val.astype(np.float64)
    # print("y_val's shape" + str(y_val.shape))

    y_train = y_train.astype(np.float64)
    y_train = y_train.to_numpy()
    #y_train = np.stack([y_train, 1 - y_train])
    #y_train = y_train.transpose()
    # print("y_train's shape" + str(y_train.shape))

    y_test = y_test.astype(np.float64)
    y_test = y_test.to_numpy()
    #y_test = np.stack([y_test, 1 - y_test])
    #y_test = y_test.transpose()
    # print("y_test's shape" + str(y_test.shape))

    X_train = X_train.to_numpy()
    X_test = X_test.to_numpy()

    # print("X_train's shape" + str(X_train.shape))
    # the real train and test dataset
    train_dataset = TensorDataset(torch.from_numpy(X_train).float(), torch.from_numpy(y_train).long())
    val_dataset = TensorDataset(torch.from_numpy(X_val).float(), torch.from_numpy(y_val).long())
    test_dataset = TensorDataset(torch.from_numpy(X_test).float(), torch.from_numpy(y_test).long())

    class NeuralNet(nn.Module):
        def __init__(self):
            super(NeuralNet, self).__init__()
            self.n = X_train.shape[1]  # number of rows

            self.fc1 = nn.Linear(self.n, int(self.n / 2))
            self.fc2 = nn.Linear(int(self.n / 2), int(self.n / 4))
            # self.fc3 = nn.Linear(int(len(X_train)/4),int(len(X_train)/8))
            self.output = nn.Linear(int(self.n / 4), 2)

            # this are defining the layers and the hyper paramters that means the
            # conditions to compare

        def forward(self, x):
            x = F.relu(self.fc1(x))
            x = F.relu(self.fc2(x))
            # x = F.relu(self.fc3(x))
            out = self.output(
                x)  # the RelU is non linearity that removes all negative values and shouldn't be used right before softmax
            return out

        def print(self):
            return self.fc1

    def loss_function(prediction, target):
        #print("debug: ",prediction.shape, target.shape)
        #print(target)
        loss = torch.nn.CrossEntropyLoss()(prediction, target)
        return loss

    train_batch_size = 32  # number of data in a training batch.
    eval_batch_size = 32  # number of data in an batch.

    train_loader = torch.utils.data.DataLoader(train_dataset, batch_size=train_batch_size, shuffle=True)
    val_loader = torch.utils.data.DataLoader(val_dataset, batch_size=eval_batch_size, shuffle=False)
    test_loader = torch.utils.data.DataLoader(test_dataset, batch_size=eval_batch_size, shuffle=False)

    def train(epoch, model, train_loader, optimizer, device):

        # activate the training mode
        model.train()

        torch.set_grad_enabled(True)

        total_loss = 0
        correct = 0

        # iteration over the mini-batches
        for batch_idx, (data, target) in enumerate(train_loader):
            # transfer the data on the chosen device
            data, target = data.to(device), target.to(device)
            target = target.float()
            target = torch.cat([torch.ones_like(target) - target, target], dim=1)
            #print("debug:", target.shape, data.shape)
            #raise ValueError

            # reinitialize the gradients to zero
            optimizer.zero_grad()

            # forward propagation on the data
            prediction = model(data)

            # compute the loss function w.r.t. the targets
            loss = loss_function(prediction, target)

            # execute the backpropagation
            loss.backward()

            # execute an optimization step
            optimizer.step()

            # accumulate the loss
            total_loss += loss.item() * len(data)
            # we multiply by the length of the data incase the last minibatch is smaller we are 'denormalizing' so that a small batch isn't overweighted

            # compute the number of correct predictions
            _, pred_classes = torch.max(prediction, dim=1)
            _, target_classes = torch.max(target, dim=1)

            # print(pred_classes)
            # print(type(pred_classes))
            # print(pred_classes.shape)

            # print(prediction.shape)
            correct += int(pred_classes.eq(target_classes).sum().item())
            # print('pred', prediction)
            # print('target', target)
            # print('correct', correct)
            # print(prediction.shape)

        # compute the average loss per epoch
        mean_loss = total_loss / len(train_loader.dataset)

        # compute the accuracy
        acc = correct / len(train_loader.dataset)

        print('Train Epoch: {}   Avg_Loss: {:.5f}   Acc: {}/{} ({:.3f}%)'.format(
            epoch, mean_loss, correct, len(train_loader.dataset),
            100. * acc))

        # return the average loss and the accuracy
        return mean_loss, acc

    # Evaluation Procedure
    def evaluate(model, eval_loader, device):

        # activate the evaluation mode
        model.eval()

        total_loss = 0
        correct = 0

        with torch.no_grad():
            # we don't need to compute the gradient graph it using too much compuatation
            # iterate over the batches
            for batch_idx, (data, target) in enumerate(eval_loader):
                # transfer the data on the chosen device
                data, target = data.to(device), target.to(device)
                target = target.float()
                target = torch.cat([torch.ones_like(target) - target, target], dim=1)
                # forward propagation on the data
                prediction = model(data)

                # compute the loss function w.r.t. the targets
                loss = loss_function(prediction, target)

                # accumulate the loss
                total_loss += loss.item() * len(data)

                # compute the number of correct predictions en sortie)
                _, pred_classes = torch.max(prediction, dim=1)
                _, target_classes = torch.max(target, dim=1)

                correct += int(pred_classes.eq(target_classes).sum().item())

                # compute the average loss per epoch
        mean_loss = total_loss / len(eval_loader.dataset)

        # compute the accuracy
        acc = correct / len(eval_loader.dataset)

        print('Eval:  Avg_Loss: {:.5f}   Acc: {}/{} ({:.3f}%)'.format(
            mean_loss, correct, len(eval_loader.dataset),
            100. * acc))

        # return the average loss and the accuracy
        return mean_loss, acc

    def save_model(epoch, model, path='./'):

        # creating the file name indexed by the epoch value
        filename = path + 'neural_network_{}.pt'.format(epoch)

        # saving the model parameters
        torch.save(model.state_dict(), filename)

        return model

    def load_model(epoch, model, path='./'):

        # creating the file name indexed by the epoch value
        filename = path + 'neural_network_{}.pt'.format(epoch)

        # loading the parameters of the saved model
        model.load_state_dict(torch.load(filename))

        return model
    # maximum number of epoch
    numEpochs = 15

    # Saving frequency
    checkpoint_freq = 10

    # Directory for data backup - save the model during training
    path = './'

    # Accumulators of average losses obtained per epoch to visualize training curve
    train_losses = []
    val_losses = []

    # Performance accumulators per epoch
    train_accuracies = []
    val_accuracies = []

    # Model definition
    neural_net = NeuralNet()

    # Load the model on the chosen device
    neural_net = neural_net.to(device)

    # Optimizer definition
    optimizer = optim.Adam(neural_net.parameters(), lr=0.001)
    # optimizer = optim.SGD(neural_net.parameters(), lr=0.001)

    # Learning loop
    for epoch in tqdm(range(1, numEpochs + 1)):

        # train the model with the train dataset
        # inner loop one step
        train_loss, train_acc = train(epoch, neural_net, train_loader, optimizer, device)

        # evaluate the model with the validation dataset
        val_loss, val_acc = evaluate(neural_net, val_loader, device)

        # Save the losses obtained
        train_losses.append(train_loss)
        val_losses.append(val_loss)

        # Save the performances
        train_accuracies.append(train_acc)
        val_accuracies.append(val_acc)

        # Checkpoint
        if epoch % checkpoint_freq == 0:
            save_model(epoch, neural_net, path)

    # Save the model at the end of the training
    save_model(numEpochs, neural_net, path)

    print("\n\n\nOptimization ended.\n")

    # Activate the evaluation mode
    neural_net = neural_net.eval()

    # Select the first 10 data points of the validation set
    data, target = test_dataset[0:10]
    data = data.to(device)

    # Executing the neural network
    output = neural_net(data)  # equivalent to neural_net.forward(data)

    # Transform the output into a probability distribution with a softmax function
    output_proba = F.softmax(output, dim=1)

    # Print the probability
    # print(output_proba)

    # For each example, retrieve the class with the highest probability.
    _, prediction = torch.max(output_proba, dim=1)

    # print("Model predictions")
    # print(prediction)

    # print("Targets")
    # print(target)

    valid_loss, valid_acc = evaluate(neural_net, val_loader, device)
    test_loss, test_acc = evaluate(neural_net, test_loader, device)
    print('valid acc', valid_acc)
    print('valid loss', valid_loss)
    print('test acc', test_acc)
    print('test loss', test_loss)

    return (train_acc, valid_acc, test_acc)