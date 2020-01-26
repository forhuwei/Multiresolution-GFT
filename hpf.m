function [value,cut] = hpf(sim_mat,source,sink);
%   Hochbaum's Pseudo-flow (HPF) Algorithm Matlab implementation          
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
%   The HPF algorithm for finding Minimum-cut in a graph is described in:                                         
%   [1] D.S. Hochbaum, "The Pseudoflow algorithm: A new algorithm for the 
%   maximum flow problem", Operations Research, 58(4):992-1009,2008.      
%                                                                         
%   The algorithm was found to be fast in theory (see the above paper)    
%   and in practice (see:                                                 
%   [2] D.S. Hochbaum and B. Chandran, "A Computational Study of the      
%   Pseudoflow and Push-relabel Algorithms for the Maximum Flow Problem,  
%   Operations Research, 57(2):358-376, 2009.                             
%  
%   and                                                                   
%                                                                         
%   [3] B. Fishbain, D.S. Hochbaum, S. Mueller, "Competitive Analysis of  
%   Minimum-Cut Maximum Flow Algorithms in Vision Problems,               
%   arXiv:1007.4531v2 [cs.CV]                                             
%                                                                         
%   Usage: Within Matlab environment:                                     
%   [value,cut] = hpf(sim_mat,source,sink);                               
%                                                                         
%   INPUTS                                                                
%                                                                   
%   sim_mat - similarity matrix - a_{i,j} is the capacity of the arc (i,j)
%             a_{i,j} are non-negatives; the self-similarities (diagonal  
%             values) are zero.                                           
%             the sim_mat should be sparse (see Matlab's help)            
%   source - The numeric label of the source node                         
%   sink   - The numeric label of the sink node                           
%                                                                         
%   OUTPUTS                                                               
%                                                                  
%   value - the capacity of the cut                                       
%   cut   - the source set (see [1]), where x_i = 1, if i \in S ; 0 o/w   
%                                                                         
%   Set-up                                                                
%                                                                   
%   Uncompress the MatlabHPF.zip file into the Matlab's working directory 
%   The zip file contains the following files:                            
%   hpf.c - source code                                                   
%   hpf.mexmaci - The compiled code for Mac OS 10.0.5 (Intel)/ Matlab     
%                 7.6.0.324 (R2008a).                                     
%   hpf.mexw32  - The compiled code for Windows 7 / Matlab 7.11.0.584     
%                 (R2010b).                                               
%   demo_general - Short Matlab code that generates small network and     
%                  computes the minimum flow                              
%   demo_vision - Short Matlab code that loads a Multiview reconstruction 
%                 vision problem (see: [3]) and computes its minimum cut. 
%   gargoyle-smal.mat - The vision problem.                               
%                                                                         
%   When using this code, please cite:                                    
%   References [1], [2] and [3] above and:                                
%   B. Fishbain and D.S. Hochbaum, "Hochbaum's Pseudo-flow Matlab         
%   implementation", http://riot.ieor.berkeley.edu/riot/Applications/     
%   Pseudoflow/maxflow.html                                               
