/*
 * Needleman-Wunsch Algorithm for Sequence Alignment
 * Author: Francesco Borando
 * Date: May 2022
 * University of Milan, Department of Physics
 *
 * This C++ code implements the Needleman-Wunsch algorithm for global sequence alignment.
 * The algorithm is used to find the optimal alignment between two sequences of symbols as defined in [1]. It employs dynamic programming to compute
 * the optimal alignment score and traceback to determine the aligned sequences.
 *
 * References:
 * - [1] Needleman, S. B., & Wunsch, C. D. (1970). A general method applicable to the search for similarities
 *   in the amino acid sequence of two proteins. Journal of Molecular Biology, 48(3), 443–453.
 *
 * Usage:
 * - Compile and run the code, providing input sequences as argv[1] and argv[2].
 * - Adjust parameter scores as desired.
 * - The output will include the optimal alignment score and the aligned sequences.
 *
 */

#include <iostream>
#include <vector>
#include <cstring>
#include <cmath>

#define GAP_PENALTY -1
#define MISMATCH_SCORE -1
#define MATCH_SCORE 1

//Useful tools
int max3(int a, int b, int c);
int match_or_mismatch(char c1, char c2);
void printmatrix(int n, int m, int* M);
int score(char c1, char c2);

int main(int argc, char* argv[])
{
    if(!argv[1] || !argv[2])
    {
        std::cerr << "Please, insert sequences to confront:" << std::endl
                <<"• Sequence1 as argv[1]" << std::endl
                <<"• Sequence2 as argv[2]" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    
    const std::string s1 = argv[1], s2 = argv[2];
    const int n = s1.length(), m = s2.length();
    
    int M[n+1][m+1];
    //STEP 1: assign first row and column
    M[0][0] = 0;
    for (int i=1;i<n+1;i++)
    {
        M[i][0] = M[i-1][0] + GAP_PENALTY;
    }
    for (int i=1;i<m+1;i++)
    {
        M[0][i] = M[0][i-1] + GAP_PENALTY;
    }
    
    //STEP 2: Needelman-Wunsch matrix
    for (int i=1;i<n+1;i++)
    {
        for (int j=1;j<m+1;j++)
        {
            M[i][j] = max3(M[i-1][j-1] + match_or_mismatch(s1[i-1], s2[j-1]),
                          M[i][j-1] + GAP_PENALTY,
                          M[i-1][j] + GAP_PENALTY);
        }
    }
    
    //DEBUG
    #ifdef DEBUG
        printmatrix(n+1, m+1, &M[0][0]);
    #endif //DEBUG

    //STEP 3: Rebuild alignments 
    std::string A_1 = "";
    std::string A_2 = "";
    int i = n, j = m;
    while (i>0 || j>0)
    {
        if (i>0 
            && j>0 
            && (M[i][j] == M[i-1][j-1] + match_or_mismatch(s1[i-1], s2[j-1])))
        {
            A_1 = s1[i-1] + A_1;
            A_2 = s2[j-1] + A_2;
            i--;
            j--;
        }

        else if (i>0 
            && (M[i][j] == M[i-1][j] + GAP_PENALTY))
        {
            A_1 = s1[i-1] + A_1;
            A_2 = '-' + A_2;
            i--;
        }

        else 
        {
            A_1 = '-' + A_1;
            A_2 = s2[j-1] + A_2;
            j--;
        }
    }
    
    std::cout << "Optimal score alignment = " << M[n][m] << std::endl;
    std::cout << "A_1 : "  << A_1 << std::endl;
    std::cout << "A_2 : "  << A_2 << std::endl;

    return 0;
}


//Functions
//Return maximum of three integers
int max3(int a, int b, int c)
{
    if (a >= b && a >= c) return a;
    else if (b >= a && b >= c) return b;
    else return c;
}

//Evaluate if diagonal outcome of Needleman-Wunsch
int match_or_mismatch(char c1, char c2)
{
    return (c1 == c2) ? MATCH_SCORE : MISMATCH_SCORE;
}


//Print the matrix
void printmatrix(int n, int m, int* M)
{
    int count =0;
    for (int i=0;i<(n*m);i++)
    {
        //style
        if (M[i]>= 0)
        {
            std::cout << " " << M[i] << "  ";
        }
        else
        {
            std::cout << M[i] << "  ";
        }
        
        //endline
        count++;
        if (count == m)
        {
            std::cout << std::endl;
            count = 0;
        }
    }
}


