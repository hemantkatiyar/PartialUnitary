# What is the code for?

Given a partial unitary matrix, the goal is to find an easy to implement full unitary.  

# ToDo

- Add an option for decomposing without free evolution ****i.e.* The target gate could be decomposed into
  single qubit gates 
- This will require renaming of sections (better solution to only use LastSection, and define Ns=0)

# Steps
- Choose an initial value for number of gates(sections)(Ns) that will be part of decomposition. First section has
  only single qubit gates(Rxy) after that every section is of form (Revo-Rxy), where Revo is evolution under
  the natural Hamiltonian.
- Define number of times the decompostion will be optimised for a particular value of number of sections(Nt)
- (1) Randomly generate the aforementioned gates for Ns sections, define a targed value for the cost function (Ft) as well as check the value of the cost function(F). 
- (2) If F>Ft, stop, otherwise
    - Set Ncurrent=0
    - (3) Increment Ncurrent by 1, calculate Gradients and change the variables accordingly.
    - Check F, if F>Ft, stop, otherwise 
        - if Ncurrent<Nt goto (2) otherwise increment Ns by 1 and goto (1) 
        - goto (3)

