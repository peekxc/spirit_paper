#import "../math_ops.typ": * 
#import "../theorems.typ": *

= Concluding Remarks and Limitations <sec:concluding-remarks>
In summary, we have introduce various spectral relaxations of the persistent rank function, which has a wealth of interesting properties, including differentiability on the positive semi-definite cone, a matrix-free representation, and natural connections to variational constructions common in machine learning applications, such as the Tikhonov regularization and the Heat kernel. 
By focusing on coefficients in $bb(R)$, we were able to exploit the connection between the inner product spaces on cochain spaces and the theory of persistence measures, providing an avenue to introduce diffentiability to an otherwise discontinuous function. 
Moreover, by focusing on the spectral characterization of the rank function, we were able to study persistence in the setting of _iterative methods_, a pursuit which lead us to the Lanczos method of Krylov expansion. 
Surprisingly, these computational techniques turned out to pair well with the output-sensitive algorithm from @chen2011output, which we can use to compute either only $Gamma$-persistence pairs or the full diagram. As cycle representatives can be obtained once a pairing has been constructed, the iterative approach may be used to compute essentially any persistence invariant. 

== Limitations <sec:limitations>

There are few limitations to the proposed work that may prevent its practical use, depending on the problem and invariant of interest. We mention the limitations to the approaches proposed that we are aware of below. 

*Function parameterization*: Choosing the matrix function $Phi$ from @eq:lowner and tuning any of its associated hyper-parameters (i.e. $tau$ and $omega$ ) can be very application-dependent. We demonstrate this in @fig:codensity_opt ---for some values of $tau > 0$, the relaxation can be either have high variance due to instability or high bias due over smoothing. In the optimization setting, we addressed these issues via an iterative shrinkage approach; in the more general setting, these choices depend on a multitude of factors, such as the relative precision the target invariant needs to achieve, the amount of compute resources available, etc.
// more sophisticated methods could also be employed, such as proximal mappings.

*Rank function instability*: Though the persistence rank invariant is considered a "stable function" in the multi-parameter setting @cerri2013betti, this stability is generally proven by reduction to a matching-type distance between the pairings defining the rank function. In particular, the Betti number function $beta_p^ast$ and its various related invariants---such as the sequence containing the Betti numbers of the homology groups at all scales (the "_Betti sequence_")---is in fact _unstable_ with respect to the $1$-Wasserstein distance (see Theorem 2.5 in @johnson2021instability). We addressed this instability by a spectral approach, but other approaches that could be used; for example, @johnson2021instability provide a Gaussian-smoothing type stabilization procedure similar to persistence images @adams2017persistence, though this does require access to the pairings. Another approach---which depends only on the rank---is the _hierarchical stabilization_ method from @gafvert2017stable, though this remains unexplored.  

*Randomized estimation:* Computing some of the persistence invariants exactly---such as the pairing itself---poses some difficult in practice, especially when the exact rank is required. For example, in the divide-and-conquer algorithm from @chen2011output, a single wrong rank value can lead to unpredictable output. 
If a randomized estimator with fixed success probability $p > 0$ is used, then one may need to be re-run $k$ times to achieve a higher success probability $(1 - (1-p)^k)$. 
Fortunately, for any constant $delta in (0,1)$, Lemma 10 of @chen2011output shows the number of such rank computations needed to compute $Gamma$-persistent pairs is not larger than:

$ rho = rho(n, ceil(1/delta)) = 4 n lr((2 ceil(1 / delta) + log n + 1)) $

Indeed, the analysis from section 6 of @chen2011output shows that if the rank estimator has success probability $p$ and $hat(p) = 1 - p$, then for the entire algorithm to achieve an arbitrary success probability $q > 0$, every rank computation needs to be repeated at most $k$ times, where: 

$ log(1 slash hat(p))^(-1) dot log (omega) <= k, quad omega = (1 - q^(1/rho))^(-1) $

Since $log(omega) = Theta(rho)$, we have that every rank computation requires $O(log rho) = O(log n + log 1 / delta)$ repeated trials. In other words, the randomized aspect of the algorithm does not significantly affect its computational complexity, though it does increase the complexity of a practical implementation.

*Estimator efficiency:* In theory, perhaps the largest limitation of approximating the rank function via the GH Monte-carlo estimator via @eq:gh_trace_estimator is its convergence rate, which is dominated by $O(epsilon.alt^(-2))$. If $epsilon.alt$ is small, the number of iterations the estimator needs to _provably_ obtain a $(1 plus.minus epsilon.alt)$-approximation can be astronomical. 
To make matters worse, estimating the numerical rank may requires shrinking $epsilon.alt$ to the order of $O (1 slash n)$. Thus its practical utility lies in applications where $epsilon.alt$ need not be too small, i.e. only a relatively coarse approximation of $tr (f (A))$ is needed. From our view, this slow convergence and high precision combination of the GH estimator @eq:gh_trace_estimator is the only major barrier preventing the use of matrix-free techniques in more practical application settings.

In practice, the efficiency of the iterative implicit trace estimator from @sec:iterative_approx depends on a variety of nuanced factors, such as the condition number of the operator, the spectral gap, etc. Fortunately, implicit trace estimation is active area of research; recent results by Meyers et al. @meyer2021hutchpp, for example, show how the estimator @eq:gh_trace_estimator can be modified to obtain a $lr((1 plus.minus epsilon.alt))$ approximation using only $O( epsilon^(-1) dot sqrt( log (eta^(-1))) + log (eta^(-1)))$ samples using deflation techniques, though this does require additional re-orthogonalization and storage costs. Even more recently, optimal trace estimators based on the _exchangeability principle_ and the _Nyström approximation_ have been shown to not only match the $O(1 slash n^2)$ variance reduction rate of Hutch++, but to also empirically achieve much faster convergence @epperly2024xtrace. 

Interestingly, using a reduction to the Gap-Hamming problem, it was shown that achieving a $(1 plus.minus epsilon.alt)$ trace approximation with probability greater than $3 slash 4$ requires at least $Omega(epsilon.alt^(-1))$ matrix-vector queries $v |-> A v$ if the input vectors $v in bb(R)^n$ are chosen _non-adaptively_ @meyer2021hutchpp (i.e. independent of the structure of $A$). This lower bound suggests it may be possible to achieve a $O(1 plus.minus epsilon.alt)$ approximation of the rank-based invariants from @prop:spectral_rank_complexity with the same space complexity and without exact arithmetic, but no algorithm is known to the author that achieves this non-adaptive bound. 

// was shown to  can be obtained using $O lr((epsilon.alt^(- 1) log lr((eta^(- 1)))))$ samples chosen non-adaptively, though no algorithm is known by the author. 

// Alternative means of variance reduction include more classical approaches, such as covariate method or quasi monte carlo methods (for an example of this with Laplacians, see (cite)). We leave the optimization of @eq:gh_trace_estimator using such techniques as future work. 



// from require only of the order of $O (m n)$ when $log lr((eta^(- 1)))$ is treated as a small constant. This is not too surprising, as it in fact matches out the exact bound from. however we consider it a limitation nonetheless.