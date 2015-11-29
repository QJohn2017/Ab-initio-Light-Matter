# Rationale #
At present, Quantum propagation codes for light-matter (i.e single atoms and molecules) interactions tend to be shared only informally between physicists, or published then abandoned. This project aims to publicly develop a reliable code while using cutting edge numerical methods with an emphasis on efficiency and accuracy.

# Notes #
* Vulkan usage is dependent on double precision support (the spec hasn't been released yet)
* This is very much a work in progress.
* Optimisation will generally not be considered (in detail) until later when profilling can be done, so as to prevent premature optimisation.
* The various parts of the library will slowly appear.
* Zero numerical library dependencies.

I am making git commits as I work on unit tests to test specific library functionality, or when I'm relatively happy with the state of a header file.