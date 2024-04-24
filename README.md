# Multi-Cell-Line-Learning
This repository contains codes and data for the MCL paper. The following are the instructions to run the codes.

## Citation
Coming soon


## Requirements

Dependencies can be found in the `environment.yml` file.
> :warning: It is recommended to run IPOPT with HSL solvers [here](https://coin-or.github.io/Ipopt/INSTALL.html).

## Usage
### Prepare dataset

1. Split datasets into the training and testing set and fill data into the template files, `Data_input_train.xlsx` and `Data_input_test.xlsx`.
   > :warning: Notation change, $\beta$ of enzymes are equivailent to `E0d` columns in the excel file.
2. Provide feasbile solutions of each steady-state time point to create warmstart points for the algorithm. The input template can be found in `\SCL_MCL\Data\IGs\CL1\solution_CL1_1.csv`. The filename (path) should follow the following template `\SCL_MCL\Data\IGs\$cell-line-name\solution_$cell-line-name_$random-number.csv`.

### Steps to run SCL/MCL optimization codes in the `SCL_MCL` folder
1. **Solve SCL for each cell line**
   ```
   python Wrapper.py -stg 1 -ti $cell-line-name
   ```
   Run the code with the specified `$cell-line-name`.
2. **Perform 5-fold cross-validation to tune lambda in MCL**

   For each $\lambda$ value, change the input values of `-l1` in the following code in the log scale, i.e. -1 for $\lambda = 0.1$.
   ```
   python Wrapper.py -kidx $kidx -l1 $l1 -np 16 -stg 2
   ```
   Run the code with `-kidx` from 1 to 5. Then, run the following code to do the model validation on the validation set.
   ```
   python Wrapper.py -kidx $kidx -l1 $l1 -np 16 -stg 3
   ```
4. **Solve MCL on the full dataset**
   
   Based on the previous stage, pick up the optimal $\lambda$ value and run 
   ```
   python Wrapper.py -l1 $l1 -np 16 -stg 4
   ```
5. **Model Validation on the testing data set**

   Get the testing error for SCL:
   ```
   python Model_Validation.py -stg 1
   ```

   Get the testing error for MCL:
   ```
   python Model_Validation.py -l1 $l1 -stg 2
   ```
### Steps to run the final cell-line models based on MCL in the `Final_Model` folder
1. **Move MCL solutions to `input\Model_inputs_MCL.csv`**
2. **Run the steady-state simulation using the specified cell line and flux-state model**
   ```
   python ModelWrapper.py -CL CL1 -fs HF -mu_target 1.0 -gln_target 2.5 -glcUB 30 -glcLB 5 -lacUB 15 -lacLB 2.0 -step_number 20 -init_step_size 0.5 -stop_step_size 1e-1 -fname SS_sol
   ```
   Read the `ModelWrapper.py` file for the details of the input arguments.
