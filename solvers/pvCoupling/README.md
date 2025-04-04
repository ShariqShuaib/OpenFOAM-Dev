This solver presents the implementation of multiple well-known pressure–velocity coupling algorithms, applied to the steady, incompressible passive scalar transport equation.

- The solvers are implemented from scratch, without relying on OpenFOAM’s built-in implicit or explicit relaxation mechanisms.
- Manual relaxation and correction are handled explicitly within the code, for better readability and understanding the basics.
- The file `createFields.H` contains all necessary declarations and initialisations required by the solver.
