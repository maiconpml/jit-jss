#include <iostream>
#include <vector>
using namespace std;

class Job{
private:
    
public:
    unsigned M; //Machine Index
    unsigned J; //Job Index
    unsigned P; //Job Cost
    unsigned S; //Job Start Time
    unsigned E; //Job End Time
    unsigned O; //Order of the job related to the oders

    void printJob(){
        cout << "m:" << M << " j:" << J << " s:" << S << " e:" << E <<endl;
    }

};

class Solution{
private:

    void printVector(vector<Job> j, unsigned sz){
        for(int i = 0; i< sz; i++){
            j[i].printJob();
        }
    }

    void printMachineSolution(){
        for(int i = 0;i < M; i++){
            printVector(machineOrder[i],J);
        }
    }

    void printJobSolution(){
        for(int i = 0;i < J; i++){
            printVector(jobOrder[i],M);
        }
    }
    bool verifyJobsAntecedence(){

        for(int i = 0; i< J; i++){
            if(!verifyjobAntedecende(jobOrder[i])) return false;
        }
        return true;
    }

    bool verifyjobAntedecende(vector<Job> job){
        for(int i = 0; i < M-1; i++){
            if(!(job[i].E <= job[i+1].S)){
                job[i].printJob();
                job[i+1].printJob();
                return false;
            } 
        }
        return true;
    }

    void swap(unsigned i1,  unsigned i2,unsigned mI){
        Job aux = machineOrder[mI][i1];
        machineOrder[mI][i1] = machineOrder[mI][i2];
        machineOrder[mI][i2] = aux;
    }

    void generateMachineOrder(vector<Job> machine, unsigned mI/*Machine index*/){
       
    }

    void generateMachinesOrder(){
        for(int m = 0; m <M; m++){
            for(int i = 0; i< J; i++){
            machineOrder[m][i] = machineJob[m][i];      
        }

        for(int i = 0; i < J; i++){
            unsigned menor = i;
            for(int j = i; j < J; j++){
                if(machineOrder[m][j].S <= machineOrder[m][menor].S){
                    menor = j;
                }
            }
            swap(menor,i,m);
        }
        }
    }

    bool verifyMachineSchedule(vector<Job> machine){
        for(int i = 0; i<J-1; i++){
            if(!(machine[i].E <= machine[i+1].S)){
                machine[i].printJob();
                machine[i+1].printJob();
                return false;
            } 
        }
        return true;
    }
    
    bool verifyMachinesSchedule(){
        generateMachinesOrder();
         for(int i = 0; i< M; i++){
            if(!verifyMachineSchedule(machineOrder[i])) return false;
        }
        return true;
    }

public:
    vector<vector<Job>> machineJob;
    vector<vector<Job>> machineOrder;
    vector<vector<Job>> jobOrder;
    unsigned M,J;
    Solution(unsigned M, unsigned J) : M(M), J(J){
        machineJob.resize(M);
        machineOrder.resize(M);
        jobOrder.resize(J);
        
        for(int i = 0; i < M; i++){
            machineJob[i].resize(J);
            machineOrder[i].resize(J);
        }
        for(int i = 0; i < J; i++){
            jobOrder[i].resize(M);

        }
        unsigned m,j;
        for(int i = 0; i < M * J ; i++){
            cin >> m;
            cin >> j;
            cin >> machineJob[m][j].P;
            cin >> machineJob[m][j].S;
            cin >> machineJob[m][j].O;
            machineJob[m][j].M = m;
            machineJob[m][j].J = j;
            machineJob[m][j].E = machineJob[m][j].P + machineJob[m][j].S;
            jobOrder[j][machineJob[m][j].O] = machineJob[m][j];
            
        }
    };

    bool isViable(){
        bool jobAntecedence = verifyJobsAntecedence();
        bool machineorder = verifyMachinesSchedule();
        if(jobAntecedence && machineorder) cout << "OK!";
        else{
            cout << "ERROR";
            if(!jobAntecedence) cout << "-JOB STARTS BEFORE ITS PREDECESSOR ENDS";
            if(!machineorder) cout << "-TWO JOBS USING THE SAME MACHINE AT SAME TIME";
            

        }
        cout<< endl;
        return jobAntecedence && machineorder;
        
    }
};
