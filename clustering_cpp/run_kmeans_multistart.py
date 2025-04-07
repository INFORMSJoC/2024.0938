from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
import pandas as pd
import numpy as np
import argparse
import os

def read_number_from_file(file_path):
    try:
        with open(file_path, 'r') as file:
            number = file.readline().strip()
            return float(number)  # Convert the read value to a float
    except FileNotFoundError:
        print(f"The file {file_path} does not exist.")
    except ValueError:
        print("The file does not contain a valid number.")
    except Exception as e:
        print(f"An error occurred: {e}")

def write_number_to_file(file_path, number):
    try:
        with open(file_path, 'w') as file:
            file.write(str(number))  # Convert the number to a string
    except Exception as e:
        print(f"An error occurred: {e}")

if '__main__' == __name__:
    parser = argparse.ArgumentParser()
    parser.add_argument("data_file", help="Full path of the dataset", type=str)
    parser.add_argument("k", help="Number of clusters", type=int)
    parser.add_argument("n_init", help="Number of times the k-means algorithm is run with different centroid seeds", type=int)
    parser.add_argument("assignment_file", help="Full path of the assignment vector", type=str)

    args = parser.parse_args()

    X = pd.read_csv(args.data_file, sep=" ")
    # print()
    print('Input data shape:', X.shape)

    if os.path.isfile(args.assignment_file[:-4] + '_MSSC.txt'):
        best_mssc = read_number_from_file(args.assignment_file[:-4] + '_MSSC.txt')
        print("\nExisting MSSC objective: " + str(best_mssc))
    else:
        best_mssc = float('inf')

    best_clusters = []

    for i in range(0, args.n_init):
        kmeans = KMeans(n_clusters=args.k, n_init=1).fit(X)
        if kmeans.inertia_ < best_mssc:
            print("Trial: " + str(i) + "\t MSSC objective: " + str(kmeans.inertia_))
            best_mssc = kmeans.inertia_
            best_clusters = kmeans.labels_
            df = pd.DataFrame(best_clusters, columns=['cluster_id'])
            df.to_csv(args.assignment_file, sep=" ", index=False, header=False)
            write_number_to_file(args.assignment_file[:-4] + '_MSSC.txt', best_mssc)