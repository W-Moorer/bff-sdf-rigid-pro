# Theory

## Contact Atlas

For a rigid rest surface patch \(S_i^0\), the contact atlas stores a chart

\[
f_i:S_i^0 \rightarrow \Omega_i \subset \mathbb R^2.
\]

The BFF adapter is a preprocessing step. It is called when a patch is created or chart topology changes. Runtime contact queries reuse stored UV coordinates and never recompute BFF per frame.

## Rigid Transform

Rest/local coordinates are named \(X0\). World coordinates are named \(xWorld\).

\[
xWorld = R X0 + t,\qquad
xLocal = R^T(xWorld-t).
\]

SDF gradients transform as

\[
nWorld = R nLocal.
\]

## Pulled-Back Signed Gap

For a source chart point \(u\), the source position is

\[
x_A(u,t)=R_A X_A^0(u)+t_A.
\]

The target SDF is queried in target local coordinates:

\[
y_B = R_B^T(x_A(u,t)-t_B).
\]

The pulled-back gap field is

\[
g_{A\rightarrow B}(u,t)=\Phi_B^0(y_B).
\]

This creates a continuous signed field over the chart for a fixed linear source mesh.

## UV Gradient

For a UV triangle, \(X0(u)\) is affine:

\[
X0(u)=\alpha X_a+\beta X_b+\gamma X_c.
\]

The chart Jacobian is

\[
J=\frac{\partial xWorld}{\partial u}=R_A\frac{\partial X0}{\partial u}.
\]

Given target normal \(nWorld\), the chart-space gradient is

\[
\nabla_u g = J^T nWorld.
\]

## Marching Triangle

Each UV triangle stores scalar values \(g_0,g_1,g_2\). The contact region is clipped by \(g\le 0\). On an edge \((i,j)\), the zero crossing is

\[
\lambda=\frac{g_i}{g_i-g_j},\qquad
u_*=(1-\lambda)u_i+\lambda u_j.
\]

The implementation handles all-positive, all-negative, one-negative, two-negative, and exact-zero vertex cases and removes duplicate zero points.

## Contact Area

The negative clipped polygon in UV is mapped back to the source triangle surface. Its 3D area is computed by fan triangulation:

\[
A=\sum_k \frac12 \lVert (p_k-p_0)\times(p_{k+1}-p_0)\rVert.
\]

The sphere-plane validation uses

\[
A=2\pi R\delta,\qquad
r=\sqrt{2R\delta-\delta^2}.
\]

The sphere-sphere equal-radius validation uses cap height \(\delta/2\):

\[
A=\pi R\delta.
\]

## Projection Refinement

SDF queries provide fast candidate filtering and closest-point hints. Atlas-Linear refinement projects a candidate point to the exact closest point on the target linear triangle mesh. The result records convergence, face id, barycentric coordinate, normal, signed gap, and fallback counters.

Atlas-HO is reserved as an extension point for PN triangles, MLS, or subdivision surfaces. The current implementation does not claim high-order accuracy.

