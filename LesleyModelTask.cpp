#include <iostream>
#include <iomanip>
#include <list>
#include <algorithm>
#include <string>
#include "gsl_poly.h"

//Matrix size
const int size = 6;

//The element of the equation or matrix
//For example, the equation x^2 + 3*x - 10 has THREE elements: Element(2, 1), Element(1, 3), Element(0, -10)
struct Element {
    //If power > 0, we have the Lambda, possibly with a coefficient. If power = 0, we have the simple number
    int power;
    double coef;

    //Are two elements the same
    bool operator == (Element& other) {
        return this->power == other.power && this->coef == other.coef;
    }

    //Add one element to another
    void operator + (Element& other) {
        this->coef += other.coef;
    }

    //Multiply one element by another
    void operator * (Element& other) {
        this->power += other.power;
        this->coef *= other.coef;
    }

    //For element sorting - the element with greater power is assumed bigger
    bool operator < (Element& other) {
        return this->power < other.power;
    }
};

//Output function for element list (polynomial)
void PrintElements(std::list<Element> elements) {
    std::list<Element>::iterator iter = elements.begin();
    std::string temp;
    for (iter; iter != elements.end(); iter++) {
        if (iter != elements.begin() && (*iter).coef >= 0)
            temp += " + ";
        temp += std::to_string((*iter).coef);
        temp += "*x^";
        temp += std::to_string((*iter).power) + " ";
    }
    std::cout << std::setw(25) << temp;
}

//Output function for Element matrix
void PrintMatrix(std::list<Element> matrix[size][size]) {
    for (int r = 0; r < size; r++) {
        for (int c = 0; c < size; c++) {
            PrintElements(matrix[r][c]);
            std::cout << "\t";
        }
        std::cout << "\n";
    }
}

//Output function for numeric matrix
void PrintMatrix(double matrix[size][size]) {
    for (int r = 0; r < size; r++) {
        for (int c = 0; c < size; c++) {
            std::cout << std::setw(10) << matrix[r][c];
            std::cout << "\t";
        }
        std::cout << "\n";
    }
}

//If we have elements like (0 * x^2), they must be removed. Only the (0 * x^0) can stay
void ClearElements(std::list<Element> &elements) {
    std::list<Element>::iterator iter = elements.begin();
    std::list<Element> temp;
    for (iter; iter != elements.end(); iter++)
        if ((*iter).coef != 0.0 || ((*iter).coef == 0.0 && (*iter).power == 0.0)) {
            temp.push_back((*iter));
        }
    elements = temp;
}

//Combine the elements with the same power. Example: x^2 + 2*x^2 - 4*x = 3*x^2 - 4*x
std::list<Element> ListCombine(std::list<Element> toCombine) {
    std::list<Element> results;
    std::list<int> powers;
    std::list<Element>::iterator secondIter = toCombine.begin();
    std::list<Element> toChange(toCombine);
    std::list<Element>::iterator iterToChange = toChange.begin();
    for (std::list<Element>::iterator iter = toCombine.begin(); iter != toCombine.end(); iter++, iterToChange++) {
        for (secondIter = toCombine.begin(); secondIter != toCombine.end(); secondIter++) {
            if ((*iter).power == (*secondIter).power && &*iter != &*secondIter) {
                (*iterToChange).operator+(*secondIter);
            }
        }
        bool contains = false;
        for (std::list<int>::iterator it = powers.begin(); it != powers.end(); it++) {
            if ((*it) == (*iter).power) {
                contains = true;
                break;
            }
        }
        if (!contains) {
            results.push_back(*iterToChange);
            powers.push_back((*iter).power);
        }
    }
    return results;
}

//Add two cells. Example: (x^2 - 1) + (x + 3) = x^2 + x + 2
std::list<Element> Add(std::list<Element> Left, std::list<Element> Right) {
    std::list<Element> results(Left);
    std::list<int> powers;
    std::list<Element>::iterator iterL = Left.begin();
    std::list<Element>::iterator iterR = Right.begin();
    for (iterL; iterL != Left.end(); iterL++)
        powers.push_back((*iterL).power);
    iterL = Left.begin();
    for (iterR; iterR != Right.end(); iterR++) {
        std::list<int>::iterator powerIt = std::find(powers.begin(), powers.end(), (*iterR).power);
        if (powerIt == powers.end()) {
            results.push_back((*iterR));
            powers.push_back((*iterR).power);
        }
        else {
            for (iterL = results.begin(); iterL != results.end(); iterL++) {
                if ((*iterR).power == (*iterL).power) {
                    (*iterL).operator+(*iterR);
                    break;
                }
            }
        }
    }
    results.sort();
    std::list<Element> finalRes(ListCombine(results));
    return finalRes;
}

//Multiply two cells. Example: (x^2 - 1) * (x + 3) = x^3 + 3*x^2 - x - 3
std::list<Element> Multiply(std::list<Element> Left, std::list<Element> Right) {
    std::list<Element> results;
    std::list<Element>::iterator iterR = Right.begin();
    std::list<Element>::iterator iterL = Left.begin();
    for (iterR; iterR != Right.end(); iterR++) {
        std::list<Element> temp(Left);
        for (iterL = temp.begin(); iterL != temp.end(); iterL++) {
            (*iterL).operator*(*iterR);
        }
        temp.sort();
        results.merge(temp);
    }
    results.sort();
    std::list<Element> finalRes(ListCombine(results));
    return finalRes;
}

//Add two matrixes
void MatrixAdd(std::list<Element> (&matrix1)[size][size], std::list<Element> (&matrix2)[size][size], 
    std::list<Element> (&result)[size][size]) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            result[i][j] = Add(matrix1[i][j], matrix2[i][j]);
        }
    }
}

//Multiply two matrixes
void MatrixMultiply(std::list<Element> (&matrix1)[size][size], std::list<Element>(&matrix2)[size][size], 
    std::list<Element>(&result)[size][size]) {
    for (int row = 0; row < size; row++) {
        for (int col = 0; col < size; col++) {
            for (int k = 0; k < size; k++) {
                result[row][col] = Add(result[row][col], Multiply(matrix1[row][k], matrix2[k][col]));
            }
        }
    }
}

//The resulting coefficients will be put here
std::list<std::list<Element>> polynomeCoefs;

std::list<Element> CalculateDiagonalSum(std::list<Element> (&matrix)[size][size]) {
    std::list<Element> results({ Element{ 0, 0.0 } });
    for (int i = 0; i < size; i++)
        results = Add(results, matrix[i][i]);
    return results;
}

//Find the coefficients of the characteristic polynomial
void FaddeevLeVerrierAlgo(std::list<Element>(&matrix1)[size][size], std::list<Element>(&matrix2)[size][size]) {
    Element polynomeCoef = { 0, 1.0 };
    polynomeCoefs.push_front(std::list<Element>{ polynomeCoef });
    for (int i = 1; i <= size; i++) {
        //M+1
        std::list<std::list<Element>>::iterator iter = polynomeCoefs.begin();
        for (int k = 0; k < size; k++) {
            matrix2[k][k] = Add(matrix2[k][k], (*iter));
        }
        //A*M
        std::list<Element> temp[size][size];
        for (int r = 0; r < size; r++) {
            for (int c = 0; c < size; c++) {
                temp[r][c] = std::list<Element>{Element{ 0, 0.0 }};
            }
        }
        MatrixMultiply(matrix1, matrix2, temp);
        for (int r = 0; r < size; r++) {
            for (int c = 0; c < size; c++) {
                ClearElements(temp[r][c]);
            }
        }
        //Coef calculation
        Element multiplier{ 0, -1.0 / (double)i };
        std::list<Element> diagonal = CalculateDiagonalSum(temp);
        std::list<Element> coef = Multiply(std::list<Element>{multiplier}, diagonal);
        polynomeCoefs.push_front(coef);
        for (int r = 0; r < size; r++) {
            for (int c = 0; c < size; c++) {
                matrix2[r][c] = temp[r][c];
            }
        }
    }
}

int main()
{
    //Form the matrixes
    std::list<Element> array[] = { std::list<Element>{Element{0, 0.72}}, std::list<Element>{Element{0, 0.0}}, std::list<Element>{Element{0, 0.0}}, std::list<Element>{Element{1, 3.6}, Element{0, -3.6}}, std::list<Element>{Element{1, 5.1}, Element{0, -5.1}}, std::list<Element>{Element{1, 7.5}}, 
        std::list<Element>{Element{0, 0.28}}, std::list<Element>{Element{0, 0.69}}, std::list<Element>{Element{0, 0.0}}, std::list<Element>{Element{0, 0.0}}, std::list<Element>{Element{0, 0.0}}, std::list<Element>{Element{0, 0.0}}, 
        std::list<Element>{Element{0, 0.0}}, std::list<Element>{Element{0, 0.31}}, std::list<Element>{Element{0, 0.75}}, std::list<Element>{Element{0, 0.0}}, std::list<Element>{Element{0, 0.0}}, std::list<Element>{Element{0, 0.0}}, 
        std::list<Element>{Element{0, 0.0}}, std::list<Element>{Element{0, 0.0}}, std::list<Element>{Element{0, 0.25}}, std::list<Element>{Element{0, 0.77}}, std::list<Element>{Element{0, 0.0}}, std::list<Element>{Element{0, 0.0}}, 
        std::list<Element>{Element{0, 0.0}}, std::list<Element>{Element{0, 0.0}}, std::list<Element>{Element{0, 0.0}}, std::list<Element>{Element{0, 0.23}}, std::list<Element>{Element{0, 0.63}}, std::list<Element>{Element{0, 0.0}},
        std::list<Element>{Element{0, 0.0}}, std::list<Element>{Element{0, 0.0}}, std::list<Element>{Element{0, 0.0}}, std::list<Element>{Element{0, 0.0}}, std::list<Element>{Element{0, 0.37}}, std::list<Element>{Element{0, 0.0}}};
    
    std::list<Element> A[size][size];
    std::list<Element> M[size][size];
    std::list<Element> result[size][size];
    int counter = 0;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            A[i][j] = std::list<Element>{array[counter]};
            M[i][j] = std::list<Element>{ Element{ 0, 0.0 } };
            result[i][j] = std::list<Element>{Element{ 0, 0.0 }};
            counter++;
        }
    }
    std::cout << "Initial matrix:\n";
    PrintMatrix(A);
    std::cout << "\n";
    //Find the polynome coefs
    FaddeevLeVerrierAlgo(A, M);

    //Add lambdas with powers to the coefs - using multiplying
    //Example: lambda2 coef = 2x -> the element will become 2x * x^2 = 2x^3
    std::list<std::list<Element>> polynome;
    int power = 0;
    for (auto i = polynomeCoefs.begin(); i != polynomeCoefs.end(); i++) {
        polynome.push_back(Multiply((*i), std::list<Element>{Element{ power, 1.0 }}));
        power++;
    }
    //Put all elements together as they are in separate lists
    std::list<Element> finalPolynome;
    for (auto i = polynome.begin(); i != polynome.end(); i++) {
        (*i).sort();
        finalPolynome.merge((*i));
        finalPolynome.sort();
    }
    //Combine the elements with the same power
    std::list<Element> totalPolynome = ListCombine(finalPolynome);
    ClearElements(totalPolynome);
    totalPolynome.sort();
    std::cout << "Characteristic Polynomial: ";
    std::cout << "\n";
    PrintElements(totalPolynome);
    std::cout << "\n\n";

    //Check if anything's wrong with the matrix
    if (totalPolynome.size() <= size) {
        std::cout << "It seems like your matrix is singular. Can't proceed with the calculation\n";
        return -1;
    }

    //Solve the polynome
    double a[size + 1];
    double z[size * 2];
    counter = 0;
    for (auto i = totalPolynome.begin(); i != totalPolynome.end(); i++) {
        a[counter] = (*i).coef;
        counter++;
    }
    gsl_poly_complex_workspace* w = gsl_poly_complex_workspace_alloc(size + 1);
    gsl_poly_complex_solve(a, size + 1, w, z);
    gsl_poly_complex_workspace_free(w);

    //Form the equation system
    //Max lambda
    double lambda = *std::max_element(z, z + size * 2);
    printf("Lambda (or x): %+.18f\n", lambda);
    //-Lambda*E matrix
    std::list<Element> LE[size][size];
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (i == j)
                LE[i][j] = std::list<Element>{Element{ 0, -lambda }};
            else
                LE[i][j] = std::list<Element>{Element{ 0, 0.0 }};
        }
    }
    //Translate the A matrix to the numerical matrix
    std::list<Element> numericA[size][size];
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            double temp = 0.0;
            for (auto iter = A[i][j].begin(); iter != A[i][j].end(); iter++) {
                temp += (*iter).coef * pow(lambda, (*iter).power);
            }
            numericA[i][j] = std::list<Element>{Element{ 0, temp }};
        }
    }
    std::cout << "\nNumeric matrix:\n";
    PrintMatrix(numericA);
    std::cout << "\n";
    //A1 = A - LE
    std::list<Element> A1[size][size];
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            A1[i][j] = std::list<Element>{Element{ 0, 0.0 }};
        }
    }
    MatrixAdd(numericA, LE, A1);
    //Translate the A1 matrix to numeric matrix
    //Prepare the vector of zeros for solving
    double X[size];
    double coefMatrixM[size][size];
    counter = 0;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            coefMatrixM[i][j] = (*A1[i][j].begin()).coef;
        }
    }
    std::cout << "\nEquation system coefficients:\n";
    PrintMatrix(coefMatrixM);
    std::cout << "\n";
    //Solve the equation system we got
    //Unfortunately, we must go bruteforce
    bool cont = false;
    while (true) {
        if (!cont) {
            for (int i = 0; i < size; i++)
                X[i] = INT_MAX;
            int temp = 100;
            int index = -1;
            std::cout << "Which variable do you want to set? Enter the number from 0 to " << size - 1 << "\n";
            //We assume the user will enter the correct number
            std::cin >> index;
            //Set the selected X to something
            X[index] = temp;
        }
        //And calculate the rest
        bool done = true;
        for (int eq = 0; eq < size; eq++) {  //We go by equations
            int foundXs = 0;
            for (int c = 0; c < size; c++) {
                if (X[c] != INT_MAX)
                    foundXs++;
            }
            if (foundXs == size) {  //If all Xs are found, no need to continue
                done = true;
                break;
            }

            done = true;
            int maxes = 0;
            for (int el = 0; el < size; el++)
                if (coefMatrixM[eq][el] != 0 && X[el] == INT_MAX)  //Are there any not-fount Xs in the current equation?
                    maxes++;
            if (maxes == 1) {  //We process the equation only if there's only one not-found
                cont = true;
                double sum = 0.0;
                bool firstFound = false;
                int indexToFind = -1;  //Which element we're to find now
                for (int el = 0; el < size; el++) {
                    if (X[el] != INT_MAX)  //Also form the sum
                        sum += X[el] * coefMatrixM[eq][el];
                    if (X[el] == INT_MAX && coefMatrixM[eq][el] != 0.0 && !firstFound) {
                        firstFound = true;
                        indexToFind = el;
                    }
                }
                X[indexToFind] = (-sum) / coefMatrixM[eq][indexToFind];
            }
        }
        for (int i = 0; i < size; i++) {  //If there are any not-found Xs left, we must continue
            if (X[i] == INT_MAX) {
                done = false;
                break;
            }
        }
        if (!cont && !done) {  //If the system is too complex and can't be solved
            std::cout << "Sorry! Coundn't solve your system! Please select a different variable to set!\n";
        }
        else if (done)  //If the system is solved
            break;
    }

    //Print the X vector out
    for (int i = 0; i < size; i++) {
        std::cout << "X" << i << " = " << X[i] << "\n";
    }

    //Calculate H parameter
    double H = (1.0 - 1.0 / lambda) * 100.0;
    std::cout << "H: " << H << "\n";

    return 0;
}