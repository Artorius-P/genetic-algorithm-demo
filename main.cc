#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

const int GENE_LENGTH = 22;
const int POP_SIZE = 50;
const int MAX_EVOLVE = 50;
const double Pc = 0.5;
const double Pm = 0.1;

double evaluate(unsigned long x_) {
    long len = (1 << GENE_LENGTH) - 1;
    double x = -1.0 +  (3.0*x_/len);
    const double PI = 3.14159265358979323846;
    double res = x * sin(10.0 * PI * x) + 1.0;
    return res;
}

unsigned long* crossover(unsigned long p1, unsigned long p2) {
    int c = rand() % (GENE_LENGTH - 1) + 1;
    unsigned long numLeft = (1 << c) - 1;   //by this we get a number like 00000111
    unsigned long numRight = (1 << GENE_LENGTH) -1 - numLeft;   //by this we get a number like 11111000
    unsigned long newP1, newP2;
    newP1 = (p1 & numRight) | (p2 & numLeft);
    newP2 = (p2 & numRight) | (p1 & numLeft);   
    static unsigned long res[] = {newP1, newP2};
    return res;
}

void crossoverPop(unsigned long* pop) {
    for (int i = 0; i < POP_SIZE; i++) {
        double p = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
        if (p < Pc) {
            int ptr = rand() % POP_SIZE;
            unsigned long* newP1P2 = crossover(pop[i], pop[ptr]);
            pop[i] = newP1P2[0];
            pop[ptr] = newP1P2[1];
        }
    }
}

unsigned long mutate(unsigned long p) {
    unsigned long np;
    int c = rand() % GENE_LENGTH;
    np = p ^ static_cast<unsigned long> (1 << c);
    return np;
}

void mutatePop(unsigned long* pop) {
    double p = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
    for (int i = 0; i < POP_SIZE; i++) {
        if (p < Pm) {
            unsigned long newP = mutate(pop[i]);
            pop[i] = newP;
        }
    }
}

unsigned long* init() {
    unsigned long max = (1 << GENE_LENGTH);
    static unsigned long pop[POP_SIZE];
    for (int i = 0; i < POP_SIZE; i++) {
        pop[i] = rand() % (max);
    }

    return pop;

}

double getSum(unsigned long* pop) {
    double sum = 0;
    for (int i = 0; i < POP_SIZE; i++) {
        sum += evaluate(pop[i]);
    }
    return sum;
}

unsigned long getBestIndv(unsigned long* pop) {
    int maxPtr = 0;
    double max = -1e10;
    for (int i = 0; i < POP_SIZE; i++) {
        double eval = evaluate(pop[i]);
        if (max < eval) {
            max = eval;
            maxPtr = i;
        }
    }
    return pop[maxPtr];
}

unsigned long* selectPop(unsigned long* pop) {
    int ca = 0;
    static unsigned long newPop[POP_SIZE];
    double fitness[POP_SIZE];
    double sumFitness = getSum(pop);
    double begin[POP_SIZE];
    double end[POP_SIZE];
    double start = 0;
    unsigned long lastBest = getBestIndv(pop);
    double lastEval = evaluate(lastBest);
    double min = 100;
    for (int i = 0; i <POP_SIZE; i++) {
        fitness[i] = evaluate(pop[i]);
        begin[i] = start;
        end[i] = fitness[i]/sumFitness + start; 
        start = end[i];
    }

    while (ca < POP_SIZE) {
        bool flag = false;
        double p = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
        for (int i = 0; i <POP_SIZE; i++) {
            if (p > begin[i] && p < end[i]) {
                newPop[ca] = pop[i];
                flag = true;
                break;
            }
        }
        if (!flag) {
            newPop[ca] = getBestIndv(pop);
        }
        ca++;
    }
    // replace the worst individual with last best individual
    for (int i = 0; i < POP_SIZE; i++) {
        double tmp = evaluate(newPop[i]);
        if (tmp <= min) {
            min = tmp;
        }
    }

    if (lastEval > min) {
        for (int i = 0; i < POP_SIZE; i++) {
            if (evaluate(newPop[i]) == min) {
                newPop[i] = lastBest;
            }    
        }
    }
    return newPop;
}

int main() {
    srand(8);
    int gen = 0;
    double sumEval = 0;
    static unsigned long* pop = init();
    double sum;
    unsigned long bestIndv;
    double bestEval=0;
    for(int i = 0; i < MAX_EVOLVE; i++) {
        pop = selectPop(pop);
        crossoverPop(pop);
        mutatePop(pop);
        unsigned long curBestIndv = getBestIndv(pop);
        if (bestEval < evaluate(curBestIndv)) {
            bestIndv = curBestIndv;
            bestEval = evaluate(curBestIndv);
        }
        cout << "times" << i << ": " << "current best value:"<< evaluate(curBestIndv) << " best gene: " << curBestIndv << endl;
    }
    cout << "summary:\nbest gene: " << hex << bestIndv << "\nbest value: " << bestEval << endl;
    
    return 0;

}
