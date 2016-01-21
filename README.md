# RePairs

**DO NOT USE THIS CODE!**

At least not if you want to learn something about your data,
as opposed to learning something about the method underlying this analysis.

This is an implementation of a method suggested by Gotelli and Ulrich
in their paper [*The empirical Bayes approach as a tool to identify
non-random species associations*][gu10].
Investigating the method described there, we found that the species pairs
which were considered significant by the algorithm depends *a lot* on the
number of bins, which is an arbitrary quantity.
So the results of the approach are to be considered **unreliable**.
We published our findings in [*Problems with bins: A critical reassessment
of Gotelli and Ulrich's Bayes approach using bird data*][gagern15].

This here is the re-implementation in [R][R] which we used to obtain
our results. In many ways it is more flexible than
[the original Pairs implementation][pairs] by [Ulrich][ulrich].
In particular, the source code is available and open for inspection.
(Prof. Ulrich made the source code for Pairs available to us via
personal communication. Inspecting it we did indeed find a major bug
which affects versions up to 1.1 but should be fixed in later releases.
But that does not make the sources public the way these here are.)

[gu10]: http://dx.doi.org/10.1007/s00442-009-1474-y
[gagern15]: http://dx.doi.org/10.1016/j.actao.2015.10.003
[R]: https://www.r-project.org/
[pairs]: http://www.keib.umk.pl/pairs/?lang=en
[ulrich]: http://www.keib.umk.pl/werner-ulrich/?lang=en

## License

The code is licensed under the [GNU General Public License
version 3](http://www.gnu.org/licenses/gpl-3.0)
or any later version at your option.

But please don't turn this into an R package without good reason,
since the authors consider the approach fundamentally flawed
and therefore not suitable for practical use.
