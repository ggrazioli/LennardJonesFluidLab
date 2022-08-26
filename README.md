# Lennard-Jones Fluid Simulation Lab
Lennard-Jones fluid simulation code in Python with a Jupyter-based lab activity for students


This Lennard-Jones (LJ) fluid simulation Python codes, the accompanying Jupyter notebook, and the lab activity were written by: 

Gianmarc Grazioli, Ph.D. <br>
http://gianmarc.com/

The lab activity includes a YouTube video, making it possible to assign this activity as a fully asynchronous lab activity. The video that accompanies this lab activity can be found here:
https://youtu.be/WgyuJYh1VaA


This lab activity is intended to be used as part of an undergraduate physical chemistry lab course, but really the only prerequisite is that the student knows how to take the derivative of a polynomial. The activity does require students to write some Python code, but the video linked above explains everything, so this can absolutely be completed by students who have never written a line of code in their life. In fact, I would argue that this activity is a great way to introduce someone interested in physical chemistry and/or chemical physics to programming in Python. 

The Jupyter notebook LennardJonesFluidLab.ipynb gives some background on LJ fluid simulations, and guides students along with the activity. It is recommended that students:
1. Get Anaconda installed on your computer
2. Start the video linked above
3. Follow along with the video to point of getting the Jupyter notebook open.
4. Once you have opened the Jupyter notebook, pause the video, read all the background information in the Jupyter notebook, then resume the video
5. Pause the video at each step, reading the instructions before listening to the video explanation.

More experienced programmers that want to jump straight into using the Lennard-Jones simulation should reference the main.py file, which will plot the LJ potential, run a short LJ fluid simulation, calculate a radial distribution function from the data, and create a .xyz file of the trajectory that can be opened in molecular dynamics data visualization software like VMD or PyMOL.  

Statistical mechanics afficionados should note that the simulations are of NVE ensembles, aka microcanonical ensembles. 

Have fun! 

Dr. G

Luca Arrigo Zammataro (@lucazammataro) has also written an excellent Lennard-Jones simulation code, and although this code did not fork off of it, his code was quite helpful in the development of portions of this code. Thank you, Luca! 


