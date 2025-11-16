#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <ctype.h>
#include <string.h>

/**
 * Defines tableau data type with the tableau's name, number of rows, number of columns,
 * row number of the row containing the coefficients of the LP problem objective function,
 * column number of the b column, which contains the right-hand sides of all LP problem constraints,
 * and 2D array of coefficients for all LP problem constraints
 */
typedef struct {
    char* name;
    size_t numRows;
    size_t numColumns;
    size_t objectiveFunctionRow;
    size_t bColumn;
    float** coefficients;
    char** variables;
} Tableau_t;

/**
 * Prints all coefficients in a tableau
 */
void printTableau(Tableau_t tableau)
{
    printf("%s\n", tableau.name);
    for (size_t var = 0; var < tableau.numColumns-1; var++) {
        printf(" %s \t", tableau.variables[var]);
    }
    printf(" b\n");
    for (size_t row = 0; row < tableau.numRows; row++) {
        printf("|");
        for (size_t col = 0; col < tableau.numColumns; col++) {
            printf("%.5f ", tableau.coefficients[row][col]);
        }
        printf("|\n");
    }
    printf("\n");
}

/**
 * Multiplies a row in a tableau by a multiplier
 */
void rowMultiplicationOperation(Tableau_t tableau, size_t rowNum, float multiplier)
{
    for (size_t col = 0; col < tableau.numColumns; col++) {
        tableau.coefficients[rowNum][col] *= multiplier;
    }
}

/**
 * Adds one tableau column multiplied by a multiplier to another tableau row
 */
void rowAdditionOperation(Tableau_t tableau, size_t rowNum, float multiplier, size_t additionRowNum)
{
    for (size_t col = 0; col < tableau.numColumns; col++) {
        tableau.coefficients[rowNum][col] += (multiplier * tableau.coefficients[additionRowNum][col]);
    }
}

/**
 * Finds the tableau column that has the most negative value in its objective function row
 */
int32_t mostNegativeColumn(Tableau_t tableau)
{
    int32_t mostNegativeCol = 0;
    while (mostNegativeCol < tableau.bColumn && tableau.coefficients[tableau.objectiveFunctionRow][mostNegativeCol] >= 0) {
        mostNegativeCol++;
    }
    if (mostNegativeCol == tableau.bColumn) {
        return -1;
    }

    for (size_t col = mostNegativeCol; col < tableau.bColumn; col++) {
        if (tableau.coefficients[tableau.objectiveFunctionRow][col] < tableau.coefficients[tableau.objectiveFunctionRow][mostNegativeCol]) {
            mostNegativeCol = col;
        }
    }
    return mostNegativeCol;
}

/**
 * Using the most negative column in a tableau, finds the smallest ratio between a coefficient in the b column and
 * a coefficient in the most negative column that is in the same row
 */
int32_t smallestRatioRow(Tableau_t tableau, int32_t mostNegativeCol)
{
    int32_t smallestRow = 0;
    while (smallestRow < tableau.objectiveFunctionRow && tableau.coefficients[smallestRow][mostNegativeCol] <= 0) {
        smallestRow++;
    }
    if (smallestRow == tableau.objectiveFunctionRow) {
        return -1;
    }

    float smallestRatio = tableau.coefficients[smallestRow][tableau.bColumn] / tableau.coefficients[smallestRow][mostNegativeCol];
    for (size_t row = smallestRow; row < tableau.objectiveFunctionRow; row++) {
        if (tableau.coefficients[row][mostNegativeCol] > 0) {
            float ratio = tableau.coefficients[row][tableau.bColumn] / tableau.coefficients[row][mostNegativeCol];
            if (ratio < smallestRatio) {
                smallestRatio = ratio;
                smallestRow = row;
            }
        }
    }
    return smallestRow;
}

/**
 * Using the most negative column and smallest ratio row, row operations are done to bring the
 * most negative column into the basis of the tableau, meaning the column is brought into the
 * basic feasible solution to the LP problem
 */
void bringColumnIntoBasis(Tableau_t tableau, int32_t mostNegativeCol, int32_t smallestRatioRow)
{
    float pivot = tableau.coefficients[smallestRatioRow][mostNegativeCol];
    if (pivot != 1) {
        float multiplier = 1 / pivot;
        rowMultiplicationOperation(tableau, smallestRatioRow, multiplier);
    }

    for (size_t row = 0; row < tableau.numRows; row++) {
        if (row != smallestRatioRow) {
            float coefficient = tableau.coefficients[row][mostNegativeCol];
            float multiplier = coefficient * -1;
            rowAdditionOperation(tableau, row, multiplier, smallestRatioRow);
        }
    }
}

/**
 * For each column in the basis of the tableau, each coefficient in the objective function row
 * is cleared to zero using row operations on the objective function row
 */
void clearOutBottomRow(Tableau_t tableau)
{
    for (size_t col = 0; col < tableau.numColumns; col++) {
        float columnSum = 0;
        for (size_t row = 0; row < tableau.objectiveFunctionRow; row++) {
            columnSum += tableau.coefficients[row][col];
        }
        if (columnSum == 1) {
            int8_t pivotRow = 0;
            for (size_t row = 0; row < tableau.objectiveFunctionRow; row++) {
                if (tableau.coefficients[row][col] == 1) {
                    pivotRow = row;
                }
            }
            float coefficient = tableau.coefficients[tableau.objectiveFunctionRow][col];
            float multiplier = coefficient * -1;
            rowAdditionOperation(tableau, tableau.objectiveFunctionRow, multiplier, pivotRow);
        }
    }
}

/**
 * Performs the main iterative loop of the simplex algorithm for finding an optimal basic feasible solution.
 * First clears out the bottom tableau row, finds most negative column and smallest ratio row, and brings
 * most negative column into the basis of the tableau. It continues until no more negative columns are found
 */
void simplexAlgorithmLoop(Tableau_t tableau)
{
    printTableau(tableau);

    clearOutBottomRow(tableau);
    printf("After clearing out:\n");
    printTableau(tableau);

    int32_t mostNegCol = mostNegativeColumn(tableau);
    int32_t smallestRow = -1;
    while (mostNegCol != -1) {
        printf("Most negative column: %d\n", mostNegCol);
        smallestRow = smallestRatioRow(tableau, mostNegCol);
        printf("Row with smallest ratio: %d\n", smallestRow);
        bringColumnIntoBasis(tableau, mostNegCol, smallestRow);
        printf("\nAfter bringing column %d into the basis:\n", mostNegCol);
        printTableau(tableau);

        mostNegCol = mostNegativeColumn(tableau);
    }
    printf("Stop!\n\n");
}

/**
 * Gets the file pointer for the cplex.txt file where the LP problem is written
 */
FILE* getLPFile(void)
{
    FILE* file = fopen("cplex.txt", "r");
    if (file == NULL) {
        printf("Error: Could not open cplex.txt.\n");
        exit(EXIT_FAILURE);
    } else {
        return file;
    }   
}

/**
 * Gets the number of rows of coefficients in the LP problem
 */
uint16_t getNumRows(void)
{
    FILE* file = getLPFile();
    uint16_t numRows = 0;
    uint8_t maxLineLength = 50;
    char line[maxLineLength];
    while (fgets(line, maxLineLength, file)) {
        numRows++;
    };

    fclose(file);
    return numRows;
}

/**
 * Gets the number of columns of coefficients in the LP problem
 */
uint16_t getNumColumns(void)
{
    FILE* file = getLPFile();
    uint16_t numCols = 0;
    char chr = '\0';
    char newline = '\n';
    do {
        chr = fgetc(file);
        if (isalpha(chr) || chr == newline) {
            numCols++;
        }
    } while (chr != newline);

    fclose(file);
    return numCols;
}

/**
 * Adds all coefficients in the LP problem file to the 2D coefficient array of an original tableau 
 */
void addOriginalTableauCoefficients(Tableau_t* tableau)
{
    FILE* file = getLPFile();
    char chr = '\0';

    tableau->coefficients = calloc(tableau->numRows, sizeof(float));
    for (size_t row = 0; row < tableau->numRows; row++) {
        tableau->coefficients[row] = calloc(tableau->numColumns, sizeof(float));
        for (size_t col = 0; col < tableau->numColumns; col++) {
            bool seenVariable = false;
            do {
                chr = fgetc(file);
                if (!isalnum(chr)) {
                    seenVariable = false;
                } else if (isalpha(chr)) {
                    seenVariable = true;
                }
            } while (!isdigit(chr) || seenVariable);

            float coefficient = (float)(chr - '0');
            tableau->coefficients[row][col] = coefficient;
        }
    }

    fclose(file);
}

/**
 * Adds all variable names found in the LP problem file to the variables array of an original tableau
 */
void addOriginalVariables(Tableau_t* tableau)
{
    FILE* file = getLPFile();
    tableau->variables = calloc(tableau->numColumns-1, sizeof(char*));
    char variable[50];
    size_t stringSize = 0;
    char chr = '\0';
    bool seenVariable = false;
    for (size_t var = 0; var < tableau->numColumns-1; var++) {
        do {
            chr = fgetc(file);
            if (isalpha(chr)) {
                seenVariable = true;
            }
            if (seenVariable && isalnum(chr)) {
                variable[stringSize++] = (char) chr;
            }
        } while (isalnum(chr));
        variable[stringSize] = '\0';
        tableau->variables[var] = calloc(stringSize, sizeof(char));
        strcpy(tableau->variables[var], variable);

        stringSize = 0;
        seenVariable = false;
        do {
            chr = fgetc(file);
        } while (!isalnum(chr));
    }

    fclose(file);
}

/**
 * Using an original tableau, adds all coefficients to the 2D coefficient array of an auxiliary tableau
 * as well as auxiliary columns that form an initial basic feasible solution to start the simplex algorithm
 */
void addAuxiliaryTableauCoefficients(Tableau_t* original, Tableau_t* auxiliary)
{
    size_t colWithAuxVariable = original->bColumn;

    auxiliary->coefficients = calloc(auxiliary->numRows, sizeof(float));
    for (size_t row = 0; row < auxiliary->objectiveFunctionRow; row++) {
        auxiliary->coefficients[row] = calloc(auxiliary->numColumns, sizeof(float));
        for (size_t col = 0; col < auxiliary->numColumns; col++) {
            float coefficient = 0;

            if (col < original->bColumn) {
                coefficient = original->coefficients[row][col];
            } else if (original->bColumn <= col && col < auxiliary->bColumn) {
                if (col == colWithAuxVariable) {
                    coefficient = 1;
                }
            } else if (col == auxiliary->bColumn) {
                coefficient = original->coefficients[row][original->bColumn];
            }

            auxiliary->coefficients[row][col] = coefficient;
        }
        colWithAuxVariable++;
    }
    
    auxiliary->coefficients[auxiliary->objectiveFunctionRow] = calloc(auxiliary->numColumns, sizeof(float));
    for (size_t col = 0; col < auxiliary->numColumns; col++) {
        float coefficient = 0;

        if (original->bColumn <= col && col < auxiliary->bColumn) {
            coefficient = 1;
        }

        auxiliary->coefficients[auxiliary->objectiveFunctionRow][col] = coefficient;
    }
}

/**
 * Using an original tableau, adds all variable names to the variables array of an auxiliary tableau as well as auxiliary variable names
 */
void addAuxiliaryVariables(Tableau_t* original, Tableau_t* auxiliary)
{
    auxiliary->variables = calloc(auxiliary->numColumns-1, sizeof(char*));
    for (size_t var = 0; var < auxiliary->numColumns-1; var++) {
        if (var < original->bColumn) {
            size_t stringSize = sizeof(original->variables[var]) / sizeof(original->variables[var][0]);
            auxiliary->variables[var] = calloc(stringSize, sizeof(char));
            strcpy(auxiliary->variables[var], original->variables[var]);
        } else {
            char variableNum = '0' + (var + 1);
            char variable[] = {'x', variableNum, '\0'};
            size_t stringSize = sizeof(variable) / sizeof(variable[0]);
            auxiliary->variables[var] = calloc(stringSize, sizeof(char));
            strcpy(auxiliary->variables[var], variable);
        }
    }
}

/**
 * Frees all allocated memory for a tableau
 */
void freeTableauCoefficients(Tableau_t tableau)
{
    for (size_t row = 0; row < tableau.numRows; row++) {
        free(tableau.coefficients[row]);
    }
    free(tableau.coefficients);

    for (size_t var = 0; var < tableau.numColumns; var++) {
        free(tableau.variables[var]);
    }
    free(tableau.variables);
}

/**
 * Creates an original tableau for the LP problem
 */
Tableau_t createOriginalTableau(void)
{
    uint16_t rows = getNumRows();
    uint16_t cols = getNumColumns();
    Tableau_t originalTableau = {
        .name = "Original Tableau",
        .numRows = rows,
        .numColumns = cols,
        .objectiveFunctionRow = rows - 1,
        .bColumn = cols - 1
    };
    addOriginalVariables(&originalTableau);
    addOriginalTableauCoefficients(&originalTableau);

    return originalTableau;
}

/**
 * Creates an auxiliary tableau for the LP problem
 */
Tableau_t createAuxiliaryTableau(Tableau_t originalTableau)
{
    uint8_t auxRows = originalTableau.numRows;
    uint8_t auxCols = originalTableau.numColumns + (originalTableau.numRows - 1);
    Tableau_t auxTableau = {
        .name = "Auxiliary Tableau",
        .numRows = auxRows,
        .numColumns = auxCols,
        .objectiveFunctionRow = auxRows - 1,
        .bColumn = auxCols - 1
    };
    addAuxiliaryVariables(&originalTableau, &auxTableau);
    addAuxiliaryTableauCoefficients(&originalTableau, &auxTableau);

    return auxTableau;
}

/**
 * After all auxiliary columns in the auxiliary tableau have been brought out of its basis,
 * the auxiliary tableau coefficients (excluding the auxiliary columns) are copied over to
 * the original tableau
 */
void reassignOriginalTableauCoefficients(Tableau_t original, Tableau_t auxiliary)
{
    for (size_t row = 0; row < original.objectiveFunctionRow; row++) {
        for (size_t col = 0; col < original.bColumn; col++) {
            original.coefficients[row][col] = auxiliary.coefficients[row][col];
        }
        original.coefficients[row][original.bColumn] = auxiliary.coefficients[row][auxiliary.bColumn];
    }
}

/**
 * After an optimal basic feasible solution is found in the original tableau,
 * the optimal solution and optimal value is printed
 */
void printSolution(Tableau_t tableau)
{
    printf("Optimal solution:\n");
    for (size_t col = 0; col < tableau.bColumn; col++) {
        float variableValue = 0;
        for (size_t row = 0; row < tableau.objectiveFunctionRow; row++) {
            if (tableau.coefficients[row][col] == 1) {
                variableValue = tableau.coefficients[row][tableau.bColumn];
            }
        }
        printf("Variable %s: %f\n", tableau.variables[col], variableValue);
    }
    printf("\n");

    float optimalValue = tableau.coefficients[tableau.objectiveFunctionRow][tableau.bColumn] * -1;
    printf("Optimal Value: %f\n", optimalValue);
}

/**
 * Main program.
 * Creates original and auxiliary tableaus, performs simplex algorithm on the auxiliary tableau,
 * basic feasible solution from the auxiliary tableau is copied over to the original tableau,
 * performs simplex algorithm on the original tableau, prints optimal basic feasible solution,
 * and frees all allocated memory
 */
int main(void)
{
    Tableau_t originalTableau = createOriginalTableau();
    Tableau_t auxTableau = createAuxiliaryTableau(originalTableau);

    simplexAlgorithmLoop(auxTableau);

    reassignOriginalTableauCoefficients(originalTableau, auxTableau);

    simplexAlgorithmLoop(originalTableau);

    printSolution(originalTableau);

    freeTableauCoefficients(auxTableau);
    freeTableauCoefficients(originalTableau);

    return EXIT_SUCCESS;
}