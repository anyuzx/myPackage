//
//  This is a lattice random walk config generator using pivot algorithm.
//  can be self-avoid or not self-avoid
//  lattice_SAW(chain array, number of steps(monomers), step length(bond), self-avoid or no-self-avoid)
//
//  Created by guangshi on 8/25/14.
//  Copyright (c) 2014 Guang Shi. All rights reserved.
//

#include <math.h>
#include <random>

using namespace std;

void c_lattice_chain(double* chain, int N, double l0, int ve, int t){
    int i, j;
    int acpt=0, duplicate;
    int pick_pivot, pick_side, pick_rotate;
    int symtry_oprtr[3][3];
    double c_x = 0, c_y = 0, c_z = 0;
    // define the rotation matrix used in pivot algorithm
    int rotate_matrix[9][3][3]={{{1,0,0},{0,0,-1},{0,1,0}},{{1,0,0},{0,-1,0},{0,0,-1}},{{1,0,0},{0,0,1},{0,-1,0}},
        {{0,0,1},{0,1,0},{-1,0,0}},{{-1,0,0},{0,1,0},{0,0,-1}},{{0,0,-1},{0,1,0},{-1,0,0}},
        {{0,-1,0},{1,0,0},{0,0,1}},{{-1,0,0},{0,-1,0},{0,0,1}},{{0,1,0},{-1,0,0},{0,0,1}}
    };
    
    // use mt19937 random engine
    mt19937 re((unsigned int)time(0));
    // set uniform integer distribution for choosing pivot point
    std::uniform_int_distribution<int> ui_pivot(1,N-2);
    // set uniform integer distribution for choosing side for rotating
    std::uniform_int_distribution<int> ui_side(1,2);
    // set unifrom interger distribution for choosing rotate matrix
    std::uniform_int_distribution<int> ui_rotate(0,8);
    
    // initialize the chain position
    // Initially, chain is set to be straight along x-axis
    for (i=0;i<N;i++){
        for (j=0;j<3;j++){
            if (j==0){
                chain[i*3+j] = double(i+1);
            }
            else{
                chain[i*3+j] = double(0);
            }
        }
    }
    
    // pivot algorithm. The number of successful rotating operations
    // need to be the same as the number of steps
    while (acpt <= t){
        // duplicate is used to check whether there are duplicate steps
        // after rotating.
        duplicate = 0;
        // randomly pick the pivot point, side and rotate matrix
        pick_pivot = ui_pivot(re);
        pick_side = ui_side(re);
        pick_rotate = ui_rotate(re);
        
        // initialize the symmetry operator matrix
        for (i=0;i<3;i++){
            for (j=0;j<3;j++){
                symtry_oprtr[i][j] = rotate_matrix[pick_rotate][i][j];
            }
        }
        
        // pick_side to rotate.
        // pick_side == 1 means that the chain right to the pivot point is for rotating
        if (pick_side == 1){
            // declare the old_chain and temp_chain
            // old_chain is kept same. temp_chain is for rotating.
            double old_chain [(pick_pivot+1)*3];
            double temp_chain [(N-pick_pivot-1)*3];
            
            // initilize old_chain and temp_chain
            for (i=0;i<N;i++){
                for (j=0;j<3;j++){
                    if (i<=pick_pivot){
                        old_chain[i*3+j]=chain[i*3+j];
                    }
                    else if(i>pick_pivot){
                        temp_chain[(i-pick_pivot-1)*3+j]=chain[i*3+j];
                    }
                }
            }
            
            // apply the rotation to the chain
            int new_chain[(N-pick_pivot-1)*3];
            for (i=0;i<N-pick_pivot-1;i++){
                for (j=0;j<3;j++){
                    new_chain[i*3+j] = (temp_chain[i*3+0]-chain[pick_pivot*3+0])*symtry_oprtr[j][0]+(temp_chain[i*3+1]-chain[pick_pivot*3+1])*symtry_oprtr[j][1]+(temp_chain[i*3+2]-chain[pick_pivot*3+2])*symtry_oprtr[j][2]+chain[pick_pivot*3+j];
                }
            }
            
            // check whether it is self-avoiding or not
            if (ve == 1){
                // looping the old_chain and new_chain to check duplicates
                for (i=0;i<N-pick_pivot-1;i++){
                    for (j=0;j<=pick_pivot;j++){
                        if(new_chain[i*3+0]==old_chain[j*3+0]&&new_chain[i*3+1]==old_chain[j*3+1]&&new_chain[i*3+2]==old_chain[j*3+2]){
                            duplicate = 1;
                            break;
                        }
                    }
                    if(duplicate == 1){break;} // there are duplicates, break the loop
                }
                if(duplicate == 0){
                    for(i=pick_pivot+1;i<N;i++){
                        for(j=0;j<3;j++){
                            chain[i*3+j] = new_chain[(i-pick_pivot-1)*3+j];
                        }
                    }
                    acpt++;
                    continue; // if there is no duplicate steps, then count the number of succesful rotating steps plus 1
                }
                else if(duplicate == 1){
                    continue; // if there are duplicate steps, then start the pivot algorithm from the beginning
                }
            }
            // no self-avoiding
            else{
                for(i=pick_pivot+1;i<N;i++){
                    for(j=0;j<3;j++){
                        chain[i*3+j] = new_chain[(i-pick_pivot-1)*3+j];
                    }
                }
                acpt++;
                continue;
            }
        }
        // pick_side == 2 means that the chain left to the pivot point is for rotating
        else if(pick_side == 2){
            double old_chain [(N-pick_pivot)*3];
            double temp_chain [(pick_pivot)*3];
            for (i=0;i<N;i++){
                for (j=0;j<3;j++){
                    if (i<pick_pivot){
                        temp_chain[i*3+j]=chain[i*3+j];
                    }
                    else if(i>=pick_pivot){
                        old_chain[(i-pick_pivot)*3+j]=chain[i*3+j];
                    }
                }
            }
            
            // apply the rotate to the chain
            int new_chain[pick_pivot*3];
            for (i=0;i<pick_pivot;i++){
                for (j=0;j<3;j++){
                    new_chain[i*3+j] = (temp_chain[i*3+0]-chain[pick_pivot*3+0])*symtry_oprtr[j][0]+(temp_chain[i*3+1]-chain[pick_pivot*3+1])*symtry_oprtr[j][1]+(temp_chain[i*3+2]-chain[pick_pivot*3+2])*symtry_oprtr[j][2]+chain[pick_pivot*3+j];
                }
            }
            
            if (ve == 1){
                for (i=0;i<pick_pivot;i++){
                    for (j=0;j<N-pick_pivot;j++){
                        if(new_chain[i*3+0]==old_chain[j*3+0]&&new_chain[i*3+1]==old_chain[j*3+1]&&new_chain[i*3+2]==old_chain[j*3+2]){
                            duplicate = 1;
                            break;
                        }
                    }
                    if(duplicate == 1){break;}
                }
                if(duplicate == 0){
                    for(i=0;i<pick_pivot;i++){
                        for(j=0;j<3;j++){
                            chain[i*3+j] = new_chain[i*3+j];
                        }
                    }
                    acpt++;
                    continue;
                }
                else if(duplicate == 1){
                    continue;
                }
            }
            else{
                for(i=0;i<pick_pivot;i++){
                    for(j=0;j<3;j++){
                        chain[i*3+j] = new_chain[i*3+j];
                    }
                }
                acpt++;
                continue;
            }
        }
    }
    
    // calculate the center of the walk(chain)
    for (i=0;i<N;i++){
        c_x = c_x + chain[i*3];
        c_y = c_y + chain[i*3+1];
        c_z = c_z + chain[i*3+2];
    }

    c_x = c_x/N; 
    c_y = c_y/N; 
    c_z = c_z/N;
    int c_r[3] = {int(c_x),int(c_y),int(c_z)};
    
    // loop over the chain. move the chain to center and mulitply by step size.
    for (i=0;i<N;i++){
        for (j=0;j<3;j++){
            chain[i*3+j] = l0*(chain[i*3+j]-c_r[j]);
        }
    }
    return;
}