# Purpose
The purpose of the Bachelor Thesis was to dealing with the theory, i.e. the mathematical basics and the programming of a simulation.
In the simulation, the associated joint angles are calculated based on a given TCP trajectory, i.e. the inverse kinematics of a robot arm with 6 degrees of freedom.

# Simulation
In the following you can see the output for a given Tool-Center-Point (TCP) trajectory which maps a Lissajouis figure.
On the left side you can see the robot and the drawn trajectory, which visualizes the course of the 3 translational coordinates of the trajectory. 
The screenshot was taken in the middle of the simulation.
On the right side you can see the course of the 3 Euler angles of the TCP which is part of the trajectory so you have in sum 6 degrees of freedom.

![image](https://github.com/user-attachments/assets/1f2cd97e-3c3d-4947-83ef-e4edff7693a4)

In the following you can see the output for the auto mode, where you can set pick up coordinates with it's euler angles and the drop of coordinates with it's euler angles.
Based on both coordinates the TCP Trajectory is calculated. The calculation was a part of the Thesis as well. 
The auto mode is supposed to simulate a real application, e.g. a robot arm in a production hall that picks up a part and places it somewhere else.

![image](https://github.com/user-attachments/assets/f41e4c85-0b96-49c1-a272-83444b60c9c0)

