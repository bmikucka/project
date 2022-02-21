
# Input: x-fold validation



# This is a BALANCED set so I'm not doing any balancing
# Note that CSV2arff will do the balancing for you by using -limit

# Weka 
export WEKA=/home/amartin/weka-3-8-3
export CLASSPATH="$WEKA/weka.jar"

# Split up the demo csv file into train and test sets
# (I know the size of the file do have hard coded these)
head -50001 demo.csv > train.csv
head -1     demo.csv > test.csv
tail -11205 demo.csv >> test.csv

# I have created 'inputs.dat' to list the fields of interest as inputs
# for the machine learning

# Convert training data to arff format
csv2arff -skip -ni inputs_updated.dat dataset train.csv >train.arff
# -skip      - skip records with missing values
# -ni        - do not convert binary inputs to nominal Boolean
# inputs.dat - file containing list of input fields
# dataset    - the name of the output field in the CSV file
# train.csv  - the input csv file
# train.arff - the output arff file
# Convert test data to arff format
csv2arff -skip -ni inputs_updated.dat dataset test.csv >test.arff

# Set parameters for the machine learning
CLASSIFIER="weka.classifiers.trees.RandomForest"
NTREE=100

# This trains, saves the trained model and tests in one go
# To train without the independent testing remove the '-T test.arff'
java $CLASSIFIER -I $NTREE -t train.arff -d train.model -T test.arff > test.out

# This tests on a pre-trained model
java $CLASSIFIER -l train.model -T test.arff > test2.out


# Cleanup
# rm train.csv test.csv train.arff test.arff