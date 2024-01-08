/*
 * Hirschberg Algorithm for Sequence Alignment
 * Author: Francesco Borando
 * Date: May 2022
 * University of Milan, Department of Physics
 *
 * This C++ code implements the Hirschberg algorithm for global sequence alignment.
 * The algorithm is an improvement over the Needleman-Wunsch algorithm in terms of space complexity.
 * It uses a divide-and-conquer strategy to achieve linear space complexity, making it more memory-efficient.
 * The Hirschberg algorithm is particularly suitable for aligning long sequences.
 *
 * References:
 * - Hirschberg, D. S. (1975). A linear space algorithm for computing maximal common subsequences.
 *   Communications of the ACM, 18(6), 341–343.
 *
 * Usage:
 * - Compile and run the code, providing input sequences as argv[1] and argv[2].
 * - Adjust parameter scores as desired.
 * - The output will include the aligned sequences.
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

//overload pair sum
std::pair<std::string, std::string> operator+(std::pair<std::string, std::string> const& one, std::pair<std::string, std::string> const& two);

//sum_vectors: sum vector function
std::vector <int> sum_vectors(const std::vector<int>& v1, const std::vector<int>& v2);

//NWScore: return last line of score matrix
std::vector<int> NWScore(const std::string& X, const std::string& Y);

//NeedlemanWunsch: returns the alignment pair with standard algorithm
std::pair < std::string, std::string > NeedlemanWunsch(const std::string& X, const std::string& Y);

//argmax_element: returns position of max element in the vector argument
int argmax_element(const std::vector<int> score);

//Hirschberg: main algorithm; returns alignments-pair space-efficiently
std::pair< std::string, std::string > Hirschberg(const std::string& X, const std::string& Y);


int main(int argc, char* argv[])
{
    if(!argv[1] || !argv[2])
    {
        std::cerr << "Please, insert sequences to confront:" << std::endl
                <<"• Sequence1 as argv[1]" << std::endl
                <<"• Sequence1 as argv[2]" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    
    const std::string s1 = argv[1], s2 = argv[2];
    const int n = s1.length(), m = s2.length();
    
    std::pair<std::string, std::string> ZWpair = Hirschberg(s1,s2);
    std::cout << ZWpair.first << std::endl << ZWpair.second << std::endl;
     
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

std::vector<int> NWScore(const std::string& X, const std::string& Y)
{
    const int n = X.length();
    const int m = Y.length();
    int Score[n+1][m+1];
    std::vector<int> Lastline;
    
    //Step 1: start from zero
    Score[0][0]=0;
    
    //Step 1.1: first row penalties
    for (int j=1;j<=m;j++)
    {
        Score[0][j] = Score[0][j-1] + GAP_PENALTY;
    }
   
    for (int i=1; i<=n;i++)
    {
        Score[1][0] = Score[0][0] + GAP_PENALTY;
        for (int j=1; j<=m;j++)
        {
            Score[1][j] = max3(
                               Score[1][j-1] + GAP_PENALTY,
                               Score[0][j] + GAP_PENALTY,
                               Score[0][j-1] + match_or_mismatch(X[i-1],Y[j-1])
                               );
        }
        
        //useless, da portare sotto!
        for (int j=0;j<=m;j++)
        {
            Score[0][j] = Score[1][j];
        }
    }
    
    for (int j=0;j<=m;j++)
    {
        Lastline.push_back( Score[1][j] );
    }
    
    return Lastline;
    
}

std::pair < std::string, std::string > NeedlemanWunsch (const std::string& X, const std::string& Y)
{
    const int n = X.length(), m = Y.length();
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
    
    //STEP 2: Needelman-Wunsch
    for (int i=1;i<n+1;i++)
    {
        for (int j=1;j<m+1;j++)
        {
            M[i][j] = max3(M[i-1][j-1] + match_or_mismatch(X[i-1], Y[j-1]),
                          M[i][j-1] + GAP_PENALTY,
                          M[i-1][j] + GAP_PENALTY);
        }
    }

    
    //STEP 3: Reconstruct alignment
    std::string A_1 = "";
    std::string A_2 = "";
    int i = n, j = m;
    while (i>0 || j>0)
    {
        if (i>0
            && j>0
            && (M[i][j] == M[i-1][j-1] + match_or_mismatch(X[i-1], Y[j-1])))
        {
            A_1 = X[i-1] + A_1;
            A_2 = Y[j-1] + A_2;
            i--;
            j--;
        }

        else if (i>0
            && (M[i][j] == M[i-1][j] + GAP_PENALTY))
        {
            A_1 = X[i-1] + A_1;
            A_2 = '-' + A_2;
            i--;
        }

        else
        {
            A_1 = '-' + A_1;
            A_2 = Y[j-1] + A_2;
            j--;
        }
    }
    
    std::pair < std::string, std::string > alignment_pair;
    alignment_pair.first = A_1;
    alignment_pair.second = A_2;
    return alignment_pair;
}


std::vector<int> sum_vectors(const std::vector<int>& v1, const std::vector<int>& v2)
{
    std::vector<int> vector_sum;
    if(v1.size() != v2.size())
    {
        std::cerr << "In vector sum: vector dimensions are not equal!" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    for (int i=0; i<v1.size();i++)
    {
        vector_sum.push_back(v1[i] + v2[i]);
    }
    
    return vector_sum;
}


int argmax_element(const std::vector<int> score)
{
    int max = score[0];
    int max_index=0;
    for (int i=1; i<score.size();i++)
    {
        if(max < score[i])
        {
            max = score[i];
            max_index = i;
        }
    }
    
    return max_index;
}

//overload pair sum
std::pair<std::string, std::string> operator+(std::pair<std::string, std::string> const& one, std::pair<std::string, std::string> const& two)
{
    std::pair<std::string, std::string> pair_sum;
    pair_sum.first = one.first + two.first;
    pair_sum.second = one.second + two.second;
    return pair_sum;
}


std::pair< std::string, std::string > Hirschberg(const std::string& X, const std::string& Y)
{
    const int n = X.length();
    const int m = Y.length();
    std::pair< std::string, std::string > ZWpair;
    
    if (n==0)
    {
        for (int i=1; i<=m; i++)
        {
            ZWpair.first += '-';
            ZWpair.second += Y[i-1];
        }
        
    }
    
    else if (m==0)
    {
        for (int i=1; i<=n; i++)
        {
            ZWpair.first += X[i-1];
            ZWpair.second += '-';
        }
    }
    
    else if (n==1 || m ==1)
    {
        ZWpair = NeedlemanWunsch(X,Y);
    }
    
    else
    {
        const int xmid = n/2; //defect truncation (.5 -> .0)
        std::string X_to_xmid="",
                    X_from_xmid="",
                    X_from_xmid_rev="",
                    Y_to_ymid="",
                    Y_from_ymid="",
                    Y_rev="";
        
        //generate x[1...xmid]
        for (int i=0;i<xmid;i++)
        {
            X_to_xmid += X[i];
        }
        
        //reverse X[xmid+1 ... n] and get normal
        for (int i=0;i<(n-xmid);i++)
        {
            X_from_xmid_rev += X[n-1-i];
            X_from_xmid += X[xmid+i];
        }
        
        //reverse Y
        for (int i=1;i<=m;i++)
        {
             Y_rev += Y[m-i];
        }
        
        std::vector<int> scoreL = NWScore(X_to_xmid,Y);
        std::vector<int> scoreR = NWScore(X_from_xmid_rev,Y_rev);
        std::vector<int> scoreR_rev;
        
        //DEBUG
        #ifdef DEBUG
            std::cout << "ScoreL : ";
            for (int i=0; i<scoreL.size();i++)
            {
                std::cout << scoreL[i] << "\t";
            }
            std::cout << std::endl;
            
            //DEBUG
            std::cout << "ScoreR : ";
            for (int i=0; i<scoreR.size();i++)
            {
                std::cout << scoreR[i] << "\t";
            }
            std::cout << std::endl;
        #endif //DEBUG
        
        //reverse ScoreR
        for (int i=1;i<=scoreR.size();i++)
        {
            scoreR_rev.push_back(scoreR[scoreR.size()-i]);
        }
        
        
        const int ymid = argmax_element(sum_vectors(scoreL, scoreR_rev));
        
        //DEBUG
        #ifdef DEBUG
            std::cout << "ymid : " << ymid << std::endl;
        #endif //DEBUG
        
        //generate Y[1...ymid]
        for (int i=0;i<ymid;i++)
        {
            Y_to_ymid += Y[i];
        }
        
        //reverse X[xmid+1 ... n]
        for (int i=ymid;i<m;i++)
        {
            Y_from_ymid += Y[i];
        }
        
        
        ZWpair = Hirschberg(X_to_xmid, Y_to_ymid) + Hirschberg(X_from_xmid, Y_from_ymid);
    }
    
    return ZWpair;
}
