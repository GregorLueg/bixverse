# *bixverse package* code style

## *Why a code style?* 

*Last update: 19.03.2025*

</br></br>

If you wish to contribute to the package, please, follow this code style. It is 
not set in stone, but is just designed to generally make the maintenance of this
package easier and avoid some pitfalls of having a number of developers working
on the same code base.

## *General idea of the bixverse code style*

The general idea is to leverage modern R language features (for example using 
the new [S7 OOP system](https://github.com/RConsortium/S7) to have a more "R 
like" feeling to objects and classes, while avoiding some of the short falls of 
S3/S4), and use the [rextendr](https://github.com/extendr/rextendr) interface 
into Rust to make key computations as fast as possible. </br>
Born out of experience, the package creators want to avoid unnecessary 
complexity. Class inheritance structures that feel smart initially but make debugging
a huge pain, 8 layers of abstractions to actually find the function that is 
causing the bug or weird convoluted data structures. The motto should be at all
times: 
</br></br> 
**Keep it Simple, Stupid (KISS)**. 
</br></br> 
A good way to think about all of this is that there is a complexity budget to 
any code base and we would like to not spend too much energy on unnecessary 
abstractions over making key functionalities blazingly fast and memory-efficient.
Please, just think about this meme in doubt:

<img src="/misc/pics/stop_abstracting.png" width="350" height="504" alt="stop abstracting">

Additionally, the idea is to be quite defensive in terms of coding style and use
asserts (more to that later) where possible to validate function inputs, leverage
early returns and warnings if data is not what it should be.

## *General rules*

1. Most important one: document the code. [roxygen2](https://roxygen2.r-lib.org) 
makes this trivial, so, please, just document your functions. Provide some 
information what parameters mean, what they expect and what is being returned.
R is dynamically typed, so reasoning without going into the source code can be
in times quite painful when the documentation of a function is bad.
2. Secondly, we use [checkmate](https://mllg.github.io/checkmate/)-based asserts
for **every** function inputs. If a function expects some parameter to be double
between 0 and 1, without the option to provide `NA`, we do not want the user to 
be able to provide wrong arguments without getting an informative(!) error by
the function. This allows for quite defensive coding and one avoids a several 
minute debug session just to realise one has provided the wrong input. Also,
brings us back to point 1. If the function parameter is not covered by checkmate
, just right an extension, see [here](https://mllg.github.io/checkmate/articles/checkmate.html#extending-checkmate).
3. If the function does something beyond 'simple' transformation, aggregation of
data, renaming, plotting, etc. go to [Rust](https://www.rust-lang.org) and use the 
[rextendr](https://github.com/extendr/rextendr) interface to make computations 
go **brrrrrr** (i.e., fast). Some libraries such as [igraph](https://r.igraph.org) 
are very fast by their nature to go low level themselves, so no need to
reinvent wheels here. Nonetheless, the speed-ups you can gain from using Rust can
be incredible. Rust functions should start with *rs_*, and ideally an R wrapper
should exist to use them. Please refer to the (yet to be written) [Why Rust](/docs/why_rust.md)
section.
4. Use [data.table](https://github.com/Rdatatable/data.table) over tibble and
data.frame. *"Yeah, but I like dplyr and the tidyverse."* We get it... But the 
speed-ups, increased memory efficacy, feature richness of data.table are just 
too big to not use. data.table also inherits all of the data.frame functionality
and most dplyr code works with it, making it easy for users to jump to tidyverse
when they want to. The point of the bixverse is to be fast, so let's stick with
data.table.
5. Be explicit and defensive in the code where possible. Simple example for the
former, if you provide parameters to a function, write the parameter name. It 
makes reasoning and debugging code so much easier. Try to use meaningful variable
names. Think about future you when writing code. Will I still understand what a 
piece of code does in 12 months? If you are doubting yourself here, maybe rethink
what you wrote. For the latter, i.e., defensive, leverage early returns, asserts
(see point 2)
6. Avoid external dependencies if not absolutely necessary. The point of the 
package **is to rewrite functions from other packages into very fast, simple Rust-accelerated code** 
and reducing the (code) bloat that affects some packages in bioinformatics and
computational biology.
7. The good old for loop vs. lapply/map question... Generally speaking, our
recommendation is using `map` via [purrr](https://purrr.tidyverse.org) (or the 
equivalent parallelised versions via [furrr](https://furrr.futureverse.org), i.e., 
`future_map` derivatives) over the apply family functions. Did you not just write 
that you want to avoid external dependencies? Yeah, but map allows to make 
explicit code which is easier to reason over. `map_lgl()` is very clear that I 
will get a logical vector back. With `unlist(lapply())` it is less
obvious what is going on. For loops in R have a very bad reputation, but this is 
usually because people grow objects in memory in the loop which is a bad 
practice indeed (it is the [second circle of hell in R](https://www.burns-stat.com/pages/Tutor/R_inferno.pdf).
8. In terms of object-oriented programming, [S7](https://github.com/RConsortium/S7)
provides a way to write very R-like OOP (the methods belong to generics). For 
user-facing key methods and workflows, we recommend using this one, as most R
users will feel very familiar with it and it allows for (familiar) piping of various 
functions. In  certain cases, you might want to write a more Java or Python-like 
type of classes (i.e., encapsulated OOP), which gives you a lot of control in
terms of public vs. private methods and attributes. We recommend [R6](https://r6.r-lib.org/articles/Introduction.html)
here, but be aware that the average R users might find the R6 classes not very
intuitive. Inheritance can be quite useful in certain cases to abstract out
common generics/methods, but try to avoid deeply layered inheritance where 
possible. This is not the most complex software we are writing here, so there
should be no need for 8+ layers of inheritance.

