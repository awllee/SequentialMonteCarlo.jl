# [Integrals approximated by SMC](@id smcintegrals)

## The measures and integrals

We define a sequence of measures $\gamma_1, \ldots, \gamma_n$ by $\gamma_1 := M_1$ and for $p \in \{2,\ldots,n\}$, \\[ \gamma_p(A) := \int_{\mathsf{X}^p} \mathbf{1}_A(x_p) M_1({\rm d}x_1) \prod_{q=2}^p G_{q-1}(x_{q-1}) M_q(x_{q-1}, {\rm d}x_q), \qquad A \in \mathcal{X}. \\] We also define a second sequence of measures $\hat{\gamma}_1,\ldots,\hat{\gamma}_n$ by $\hat{\gamma}_p := \gamma_p \cdot G_p$ for $p \in \{1,\ldots,n\}$.

Letting $1$ denote the constant function $x \mapsto 1$, the normalizing constants associated with each $\gamma_p$ and $\hat{\gamma}_p$ are $\gamma_p(1)$ and $\hat{\gamma}_p(1)$, respectively. We define \\[ Z_p := \gamma_p(1), \qquad \hat{Z}_p := \hat{\gamma}_p(1), \qquad p\in \\{1,\ldots,n\\}. \\] This allows us to define their normalized, probability measure counterparts \\[ \eta_p := \gamma_p / Z_p, \qquad \hat{\eta}_p := \hat{\gamma}_p / \hat{Z}_p, \qquad p \in \\{1,\ldots,n\\}. \\]

We are now able to state which integrals SMC algorithms are built to approximate: these are simply $\gamma_p(f)$, $\hat{\gamma}_p(f)$, $\eta_p(f)$ or $\hat{\eta}_p(f)$ for some $p \in \{1,\ldots,n\}$ and some function $f:\mathsf{X}\rightarrow \mathbb{R}$.

We note from their definitions that $Z_1 = 1$ and $Z_p = \hat{Z}_{p-1}$ for $p \in \{2,\ldots,n\}$. It follows that any of the above measures can be written in terms of the sequences $\hat{Z}_1, \ldots, \hat{Z}_n$, $\eta_1, \ldots, \eta_n$, and $\hat{\eta}_1, \ldots, \hat{\eta}_n$. For example, $\gamma_p = \hat{Z}_{p-1} \eta_p$ for $p \in \{2, \ldots, n\}$. Therefore, it is sufficient to define SMC approximations of $\hat{Z}_p$, $\eta_p$ and $\hat{\eta}_p$ for $p \in \{1, \ldots, n \}$.

## [Recursive definition](@id recursivedefmeasures)

We can alternatively define the sequence of measures mentioned above recursively by $\eta_1=\gamma_1:=M_1$, \\[ \eta_p \propto \gamma_p := (\gamma_{p-1} \cdot G_{p-1}) M_p, \qquad p \in \\{2, \ldots, n\\}, \\] and for $p\in \{1,\ldots,n\}$, $\hat{\eta}_p \propto \hat{\gamma}_p := \gamma_p \cdot G_p$. It follows that one can write $\hat{Z}_1 = \eta_1(G_1)$ and $\hat{Z}_p = \hat{Z}_{p-1} \eta_p(G_p)$ for $p \in \{2, \ldots, n\}$.

We observe that $\gamma_p=\hat{\gamma}_{p-1}M_p$, and this suggests the following diagram \\[ M_1 = \gamma_1 \overset{G_1}{\longrightarrow} \hat{\gamma}_1 \overset{M_{2}}{\longrightarrow} \gamma_{2} \longrightarrow \cdots \longrightarrow \hat{\gamma}_{p-1} \overset{M_p}{\longrightarrow} \gamma_p \overset{G_p}{\longrightarrow} \hat{\gamma}_p \longrightarrow \cdots \longrightarrow \gamma_n \overset{G_n}{\longrightarrow} \hat{\gamma}_n, \\] where measures are obtained by weighting or mutation. Alternatively, we can view the sequence of probability measures using the same diagram, i.e. \\[ M_1 = \eta_1 \overset{G_1}{\longrightarrow} \hat{\eta}_1 \overset{M_{2}}{\longrightarrow} \eta_{2} \longrightarrow \cdots \longrightarrow \hat{\eta}_{p-1} \overset{M_p}{\longrightarrow} \eta_p \overset{G_p}{\longrightarrow} \hat{\eta}_p \longrightarrow \cdots \longrightarrow \eta_n \overset{G_n}{\longrightarrow} \hat{\eta}_n, \\] but where the weighting steps also involve renormalization: $\eta_p \overset{G_p}{\longrightarrow} \hat{\eta}_p$ means $\hat{\eta}_p = (\eta_p \cdot G_p) / \eta_p(G_p)$.

The measures $\hat{\gamma}_p$ and $\hat{\eta}_p$ are often referred to as the "updated" versions of $\gamma_p$ and $\eta_p$, respectively.
