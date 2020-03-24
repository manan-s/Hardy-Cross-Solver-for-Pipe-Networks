//PIPE NETWORKS SOLVER USING HARDY-CROSS TECHNIQUE
//CONVENTION: Discharge towards a node +ve, away -ve
//            In a loop, clockwise head is positive, else negative
//Nodes naming is 0,1,2,3,.....
#include <iostream>
#include <vector>
#include <cmath>
#include <list>
#include <queue>
#include <iterator>
#define g 9.80665
#define pi 3.1415926535
#define kinematic 0.00000113
//using namespace std;

double frictionFactor(double ReynoldsNumber, double RelativeRoughness,double f1){
    if(ReynoldsNumber<2500) return 64/ReynoldsNumber;
    
    double f=1/std::pow((-2*std::log10((RelativeRoughness/3.7)+(2.51/(ReynoldsNumber*std::sqrt(f1))))),2);  //Colebrook's Equation
    if (std::abs(f1-f)<0.0000001){
        return f;
    }
    else{
        return frictionFactor(ReynoldsNumber, RelativeRoughness, f);
    }
}

struct external_discharge{
    int node;
    double ext_discharge; //+ve means towards the node
};

double findExtDischarge(int i,external_discharge ext_discharge_matrix[], int len){
    for(int j=0;j<len;j++)
        if (ext_discharge_matrix[j].node==i) return ext_discharge_matrix[j].ext_discharge;
}

bool is_ext(int node, external_discharge ext_discharge_matrix[], int len){
    for(int i=0;i<len;i++)
        if (ext_discharge_matrix[i].node==node) return 1;
    return 0;
}

struct pipe{
    int n1;
    int n2;
    int from1to2;
    double length;
    double diameter;
    double epsilon;
    double discharge;
    pipe(){
        discharge=0;
        from1to2=1;
    }
};

int main(){
    int nodes,edges,loops;
    std::cin>>nodes>>edges>>loops;
    std::cout<<"Network has "<<nodes<<" nodes, "<<edges<<" pipes and "<<loops<<" loops.\n";

    std::list <int>*adjacency;
    std::list <int>*directed;
    adjacency=new std::list<int>[nodes];
    directed=new std::list<int>[nodes];

    //network pipe_network(nodes,edges,loops);
    pipe pipe_matrix[nodes][nodes];
    int to, from;

    //adding pipes and their data and filling adjacency matrix
    for(int i=0;i<edges;i++){
        
        pipe sample;

        std::cin>>sample.n1>>sample.n2;   //the from and to nodes
        std::cin>>sample.length>>sample.diameter>>sample.epsilon;

        pipe_matrix[sample.n1][sample.n2]=sample;
        sample.from1to2=-1;
        pipe_matrix[sample.n2][sample.n1]=sample;
        adjacency[sample.n1].push_back(sample.n2);
        adjacency[sample.n2].push_back(sample.n1);
        directed[sample.n1].push_back(sample.n2);
    }

    //adding loops
    std::list<int>looplist[loops];
    std::list<int>single_loop;
    for(int i=0;i<loops;i++){
        int loopsize;
        std::cin>>loopsize;
        std::cout<<"So the loopsize is: "<<loopsize<<std::endl;
        int loopnode;
        int first;
        for(int j=0;j<loopsize;j++){
            std::cin>>loopnode;
            if(j==0) first=loopnode;
            single_loop.push_back(loopnode);
            std::cout<<loopnode<<" ";
        }
        std::cout<<first<<std::endl;
        single_loop.push_back(first);
        looplist[i]=single_loop;
        single_loop.erase(single_loop.begin(),single_loop.end());
    }

    //adding external discharges
    int ext_discharge_count;
    std::cin>>ext_discharge_count;
    std::cout<<"So there are "<<ext_discharge_count<<" external discharges."<<std::endl;
    external_discharge ext_discharge_matrix[ext_discharge_count];
    
    for(int m=0;m<ext_discharge_count;m++){
        std::cin>>ext_discharge_matrix[m].node>>ext_discharge_matrix[m].ext_discharge;
        std::cout<<"Discharge at node "<<ext_discharge_matrix[m].node<<" is "<<ext_discharge_matrix[m].ext_discharge<<std::endl;
    }

    //distributing discharges
    bool discovered[nodes];
    bool visited[nodes];
    int degree[nodes];
    double temp_extq[nodes];
    //initialize
    for(int i=0;i<nodes;i++){
        discovered[i]=0;
        visited[i]=0;
        degree[i]=adjacency[i].size();
        temp_extq[i]=0;
    }

    std::queue<int>Q;
    Q.push(ext_discharge_matrix[0].node); //starting from first node in external discharge matrix
    discovered[ext_discharge_matrix[0].node]=1;
    temp_extq[ext_discharge_matrix[0].node]=ext_discharge_matrix[0].ext_discharge;
    
    while(!Q.empty()){
        std::cout<<"At node "<<Q.front()<<" having degree "<<degree[Q.front()]<<std::endl;
        int init=Q.front();
        visited[init]=1;
        discovered[init]=1;
        Q.pop();
        std::list<int> temp = adjacency[init];
        std::list<int>::iterator i;
        if(degree[init]==0) continue;
        double distributing=temp_extq[init]/degree[init];
        std::cout<<"Discharge to be distributed is "<<distributing<<std::endl;
        for(i = temp.begin(); i != temp.end(); ++i){
            
            if(discovered[*i]==0||visited[*i]==0){
                discovered[*i]=1;
                degree[*i]--;
                if(visited[*i]==0&&discovered[*i]==1) temp_extq[*i]+=(distributing+findExtDischarge(*i,ext_discharge_matrix,ext_discharge_count));
                else temp_extq[*i]+=distributing;
                Q.push(*i);
                pipe_matrix[init][*i].discharge=distributing;
                pipe_matrix[init][*i].from1to2=1;
                std::cout<<"So pipe from "<<init<<" to "<<*i<<" has discharge "<<distributing<<std::endl;
                std::cout<<"Net Discharge at "<<*i<<" is "<<temp_extq[*i];
                pipe_matrix[*i][init].discharge=-1*distributing;
                pipe_matrix[*i][init].from1to2=-1;
            }
        } 
        std::cout<<std::endl;
    }

    //Hardy-Cross method
    int numiters=10 ;
    for(int k=0;k<numiters;k++){

        std::list< std::list< int > >::iterator a;
        std::list< int >::iterator b;
        std::list< int >::iterator temp;
        for(int j=0;j<loops;j++){

            double loop_head=0;
            double loop_head_Q=0;
            int starting_node;

            for(b=looplist[j].begin(),temp=++b,b--;b!=looplist[j].end();++b){
                double roughness=pipe_matrix[*(b)][*temp].epsilon/pipe_matrix[*b][*temp].diameter;
                
                if(pipe_matrix[*b][*temp].discharge==0){
                    loop_head+=0;
                    loop_head_Q+=0;
                }
                else{
                    double reynolds=4*pipe_matrix[*b][*temp].discharge*pipe_matrix[*b][*temp].from1to2/(pi*pipe_matrix[*b][*temp].diameter*kinematic);
                    double K=(8*frictionFactor(reynolds,roughness,0.001)*pipe_matrix[*b][*temp].length)/(g*pi*pi*pipe_matrix[*b][*temp].diameter*pipe_matrix[*b][*temp].diameter*pipe_matrix[*b][*temp].diameter*pipe_matrix[*b][*temp].diameter*pipe_matrix[*b][*temp].diameter);
                    double head=K*pipe_matrix[*b][*temp].discharge*pipe_matrix[*b][*temp].from1to2*pipe_matrix[*b][*temp].discharge*pipe_matrix[*b][*temp].from1to2;
                    double err=head/(pipe_matrix[*b][*temp].discharge*pipe_matrix[*b][*temp].from1to2);
                    loop_head+=head;
                    loop_head_Q+=err;
                }
            }
            double correction=-1*loop_head/(2*loop_head_Q);
            for(b=looplist[j].begin(),temp=++b,b--;b!=looplist[j].end();++b){
                pipe_matrix[*b][*temp].discharge+=correction;
                pipe_matrix[*temp][*b].discharge-=correction;
            }
        }
    }
   for(int i=0;i<nodes;i++){
       std::list<int>temporary = directed[i];
       std::list<int>::iterator iter;
       for(iter = temporary.begin(); iter!= temporary.end(); ++i){
           std::cout<<"Discharge in Pipe "<<i<<" --> "<<*iter<<" is "<<pipe_matrix[i][*iter].discharge<<std::endl;
       }
   }
}