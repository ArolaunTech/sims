#include <iostream>
#include <vector>
#include <cmath>

#include "random.h"
#include "debug.h"

#ifndef NN_H
#define NN_H

//Utils
double clamp(double x, double min, double max) {
	if (x < min) return min;
	if (x > max) return max;
	return x;
}

//Activation functions
enum Activation {
	ACTIVATION_RELU,
	ACTIVATION_LINEAR,
	ACTIVATION_LOGISTIC,
	ACTIVATION_EXP
};

double activation_function_eval(double x, Activation type) {
	switch (type) {
		case Activation::ACTIVATION_RELU:
			if (x < 0) return 0;
			return x;
		case Activation::ACTIVATION_LINEAR:
			return x;
		case Activation::ACTIVATION_LOGISTIC:
			return 1 / (1 + std::exp(-x));
		case Activation::ACTIVATION_EXP:
			return std::exp(x);
	}

	return x;
}

double activation_function_grad(double x, Activation type) {
	switch (type) {
		case Activation::ACTIVATION_RELU:
			if (x < 0) return 0.1; //Fixes vanishing gradient problem
			return 1;
		case Activation::ACTIVATION_LINEAR:
			return 1;
		case Activation::ACTIVATION_LOGISTIC:
			return std::exp(x) / std::pow(1 + std::exp(x), 2);
		case Activation::ACTIVATION_EXP:
			return std::exp(x);
	}

	return 1;
}

//NN
struct Layer {
	std::vector<std::vector<double> > weights;
	std::vector<double> biases;
	Activation layer_activation_function;

	std::vector<double> apply_activation_function(std::vector<double> input) const {
		std::vector<double> out;
		std::size_t numValues = input.size();

		for (std::size_t i = 0; i < numValues; i++) {
			out.push_back(activation_function_eval(input[i], layer_activation_function));
		}

		return out;
	}

	std::vector<double> evaluate(std::vector<double> input, bool skip_activation = false) const {
		std::vector<double> out;
		std::size_t num_nodes_curr = weights.size();
		std::size_t num_nodes_prev = input.size();

		for (std::size_t i = 0; i < num_nodes_curr; i++) {
			double node_val = biases[i];

			for (std::size_t j = 0; j < num_nodes_prev; j++) {
				node_val += weights[i][j] * input[j];
			}

			out.push_back(node_val);
		}

		if (!skip_activation) return apply_activation_function(out);

		return out;
	}

	void mutate() {
		std::size_t num_nodes_curr = biases.size();
		std::size_t num_nodes_prev = weights[0].size();

		for (std::size_t i = 0; i < num_nodes_curr; i++) {
			biases[i] += 0.01 * normal(generator);

			for (std::size_t j = 0; j < num_nodes_prev; j++) {
				weights[i][j] += 0.01 * normal(generator);
			}
		}
	}
};

struct NN {
	std::vector<Layer> layers; //Does not include the input layer

	void log() const {
		std::size_t numLayers = layers.size();

		std::cout << "NN {\n\t- Input (" << layers[0].weights[0].size() << " neuron(s))\n";

		for (std::size_t i = 0; i < numLayers; i++) {
			std::cout << "\t- Layer (" << layers[i].biases.size() << " neuron(s)) {\n";
			std::cout << "\t\t- Activation Function: ";
			switch (layers[i].layer_activation_function) {
				case Activation::ACTIVATION_LINEAR:
					std::cout << "Linear\n";
					break;
				case Activation::ACTIVATION_RELU:
					std::cout << "RELU\n";
					break;
				case Activation::ACTIVATION_LOGISTIC:
					std::cout << "Logistic\n";
					break;
				case Activation::ACTIVATION_EXP:
					std::cout << "Exponential\n";
					break;
				default:
					std::cout << "\n";
			}

			std::size_t numNodes = layers[i].biases.size();
			for (std::size_t j = 0; j < numNodes; j++) {
				std::cout << "\t\t- Node " << j << " {\n\t\t\t- Bias: " << layers[i].biases[j] << "\n\t\t\t- Weights: {\n\t\t\t\t";

				std::size_t numNodesPrev = layers[i].weights[0].size();
				for (std::size_t k = 0; k < numNodesPrev; k++) {
					std::cout << k << ": " << layers[i].weights[j][k] << "\t";

					if (k % 5 == 4) {
						std::cout << "\n\t\t\t\t";
					}
				}

				std::cout << "\n\t\t\t}\n\t\t}\n";
			}

			std::cout << "\t}\n";
		}

		std::cout << "}\n";
	}

	std::vector<double> evaluate(std::vector<double> input) const {
		std::vector<double> out = input;

		for (Layer const & iterator : layers) {
			out = iterator.evaluate(out);
		}

		return out;
	}

	void mutate() {
		std::size_t num_layers = layers.size();

		for (std::size_t i = 0; i < num_layers; i++) {
			layers[i].mutate();
		}
	}

	void backpropagate_gradient(std::vector<double> input, std::vector<double> lastLayerGradients, double learningRate = 1e-3, double clampMagnitude = 0.1) {
		std::size_t numLayers = layers.size();

		//logDoubleVector(input);
		//std::cout << "\n";
		//logDoubleVector(output);
		//std::cout << "\n\n";

		// Forward pass
		std::vector<std::vector<double> > layerOutputs;
		std::vector<double> curr = input;

		for (Layer const & iterator : layers) {
			curr = iterator.evaluate(curr, true);

			layerOutputs.push_back(curr);

			curr = iterator.apply_activation_function(curr);
		}

		//for (std::size_t i = 0; i < numLayers; i++) {
		//	logDoubleVector(layerOutputs[i]);
		//	std::cout << "\n";
		//}
		//std::cout << "\n";

		//Backward pass
		std::vector<std::vector<double> > nodeValueGradients(numLayers);
		//std::vector<double> lastLayerGradients;
		//std::size_t lastLayerNodes = layers[numLayers - 1].biases.size();

		//for (std::size_t i = 0; i < lastLayerNodes; i++) {
		//	lastLayerGradients.push_back(output[i] - curr[i]);
		//}

		nodeValueGradients[numLayers - 1] = lastLayerGradients;

		for (int i = (int)numLayers - 2; i >= 0; i--) {
			std::vector<double> layerGradients;

			std::size_t numNodes = layers[i].biases.size();
			std::size_t numNodesAhead = layers[i + 1].biases.size();

			for (std::size_t j = 0; j < numNodes; j++) {
				layerGradients.push_back(0);

				for (std::size_t k = 0; k < numNodesAhead; k++) {
					layerGradients[j] += 
						layers[i + 1].weights[k][j] *  
						activation_function_grad(layerOutputs[i + 1][k], layers[i + 1].layer_activation_function) * 
						nodeValueGradients[i + 1][k];
				}
			}

			nodeValueGradients[i] = layerGradients;
		}

		//for (std::size_t i = 0; i < numLayers; i++) {
		//	logDoubleVector(nodeValueGradients[i]);
		//	std::cout << "\n";
		//}

		//std::cout << "\n";

		//Updates
		std::size_t numNodes = layers[0].biases.size();
		std::size_t numNodesPrev = input.size();
		for (std::size_t i = 0; i < numNodes; i++) {
			layers[0].biases[i] += clamp(
				learningRate *
				nodeValueGradients[0][i] * 
				activation_function_grad(layerOutputs[0][i], layers[0].layer_activation_function),
				-clampMagnitude, clampMagnitude
			);

			for (std::size_t j = 0; j < numNodesPrev; j++) {
				layers[0].weights[i][j] += clamp(
					learningRate * 
					input[j] *
					nodeValueGradients[0][i] * 
					activation_function_grad(layerOutputs[0][i], layers[0].layer_activation_function),
					-clampMagnitude, clampMagnitude
				);
			}
		}

		for (std::size_t i = 1; i < numLayers; i++) {
			numNodes = layers[i].biases.size();
			numNodesPrev = layers[i - 1].biases.size();

			for (std::size_t j = 0; j < numNodes; j++) {
				layers[i].biases[j] += clamp(
					learningRate *
					nodeValueGradients[i][j] * 
					activation_function_grad(layerOutputs[i][j], layers[i].layer_activation_function),
					-clampMagnitude, clampMagnitude
				);

				for (std::size_t k = 0; k < numNodesPrev; k++) {
					layers[i].weights[j][k] += clamp(
						learningRate * 
						layerOutputs[i - 1][k] *
						nodeValueGradients[i][j] * 
						activation_function_grad(layerOutputs[i][i], layers[i].layer_activation_function),
						-clampMagnitude, clampMagnitude
					);
				}
			}
		}
	}

	void backpropagate(std::vector<double> input, std::vector<double> output, double learningRate = 1e-3, double clampMagnitude = 0.1) {
		std::vector<double> actual = evaluate(input);

		std::vector<double> gradients;
		std::size_t outputNodes = output.size();
		for (std::size_t i = 0; i < outputNodes; i++) {
			gradients.push_back(output[i] - actual[i]);
		}

		backpropagate_gradient(input, gradients, learningRate, clampMagnitude);
	}
};

NN get_random_nn(std::vector<int> layer_node_counts, std::vector<Activation> layer_activations) { //Includes the input layer
	NN out;
	std::size_t num_layers = layer_node_counts.size();

	for (std::size_t i = 1; i < num_layers; i++) {
		Layer new_layer;
		new_layer.layer_activation_function = layer_activations[i];

		for (int j = 0; j < layer_node_counts[i]; j++) {
			new_layer.biases.push_back(1 * normal(generator));

			std::vector<double> node_weights;
			for (int k = 0; k < layer_node_counts[i - 1]; k++) {
				node_weights.push_back(0.1 * normal(generator));
			}

			new_layer.weights.push_back(node_weights);
		}

		out.layers.push_back(new_layer);
	}

	return out;
}

NN get_random_nn(std::vector<int> layer_node_counts, Activation layer_activations) {
	std::vector<Activation> new_layer_activations;
	std::size_t num_layers = layer_node_counts.size();

	for (std::size_t i = 0; i < num_layers - 1; i++) {
		new_layer_activations.push_back(layer_activations);
	}

	new_layer_activations.push_back(Activation::ACTIVATION_LINEAR);

	return get_random_nn(layer_node_counts, new_layer_activations);
}

bool operator<(NN const & lhs, NN const& rhs) { //Makes std::sort work
	return lhs.layers.size() < rhs.layers.size();
}

#endif