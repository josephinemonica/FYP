I copied IK.py and edited it in IK_edit.py
the Jacobian for ee w.r.t base is flipped.
It was : 
[e]
[exr]
I changed to:
[exr]
[e]
so that J q_dot=    [velocity]
                    [omega]
so that it is consistent when being concatenated to build the TOTAL jacobian

TODO
NOTE: MUST change the Jacobian of the rest of the code

-----------------------------------------------------------------------------
IDEA:
We have 2 psi (of base and of EE) careful!!
for orientation->

pyOpt for optimization in python
I have installed the package!
http://www.pyopt.org/tutorial.html

QPOEASES for finding redundancy parameters that will optimize
https://projects.coin-or.org/qpOASES

problems with joint limit

=== if theta starts from here, they will be just constant, and will not go to the middle

---

---

===
--------------------------------------------------------------------------------
CURRENT STATE:

IK_edit.py -> working code
kinematics.py -> tidy up from IK_edit.py
kinematics_clean.py -> version of kinematics.py that is PURELY CLEAN, to check the duration needed for the algorithms
IK_optimize_time.py -> optimize time for algorithms: in IK, each time we compute Jacobian and FK, but there are overlapping parts in Jacobian computation and FK computation

kinematics_JL.py -> with JOINT LIMITS

---------------------------------------------------------------------------------------------------------------------------------------
N=50 already gives good result, run in 3.49696397781 seconds
python kinematics_clean.py 
('q computed by IK algo is: ', Matrix([
[  4.66172488405933],
[   1.9897649765862],
[ 0.957053786734947],
[0.0934715234366335],
[0.0946384174341276],
[-0.694886560732354],
[  2.22747863151291],
[0.0595685724376596]]))
('the q computed would result to pose: ', Matrix([
[5.00015272966222],
[2.00006337858201],
[1.00004256932294],
[1.00000798396432],
[ 1.5000414468039],
[1.00003651653218]]))
('the error in pose is: ', Matrix([
[-0.000152729662217688],
[ -6.33785820052246e-5],
[ -4.25693229400004e-5],
[ -7.98396431833659e-6],
[ -4.14468039016125e-5],
[ -3.65165321796201e-5]]))
#========================================================================================================================
WITH JL and manipulability
q_1=Matrix([[0],[0],[2.5],[0],[0],[0+pi],[0-pi/2],[0]])   #initial posture
pose_goal=Matrix([[2],[2],[5],[1],[2],[2.1]])   #pose goal

JL const=5
manip_const=0 -> m=2.
manip_const=5 -> m=2.

JL_const=1000
manip_const=500 -> m=2.

JL_const=1000
manip_const=5000 ->m=2.

JL_const=1000
manip_const=50000 ->m=1. FAILL

