# *bixverse package* code style

The general idea is to leverage modern R language features (for example using 
the new [S7 OOP system](https://github.com/RConsortium/S7) to have a more "R 
like" feeling to objects and classes, while avoiding some of the short falls of 
S3/S4), and use the [rextendr](https://github.com/extendr/rextendr) into Rust 
to make key computations as fast as possible. </br>
Born out of experience, we would like to avoid unnecessary complexity: 
inheritance structures that feel smart but make debugging a huge pain, 8 layers
of abstractions to actually find the function that is bugging or weird convoluted
data structures. The motto should be at all times: </br> 
**Keep it Simple, Stupid (KISS)**. </br>
A good way to think about all of this is that there is a complexity budget to 
any code base and we would like to not spend too much energy on unnecessary 
abstractions over making key functionalities blazingly fast and memory-efficient.
Please, just think about this meme in doubt:

<img src="/misc/pics/stop_abstracting.png" width="350" height="504" alt="stop abstracting">

Additionally, the idea is to be quite defensive in terms of coding style and use
asserts (more to that later) where possible to validate function inputs.

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
brings us back to point 1.
3. Use [data.table](https://github.com/Rdatatable/data.table) over tibble and
data.frame. Yeah, but I like dplyr and the tidyverse. We get it... But the 
speed-ups, increased memory efficacy, feature richness of data.table are just 
too big to not use. 
4. If the function does something beyond 'simple' transformation, aggregation of
data, renaming, plotting go to [Rust](https://www.rust-lang.org) and use the 
[rextendr](https://github.com/extendr/rextendr) interface to make computations 
go *brrrrrr*. Some libraries such as [igraph](https://r.igraph.org) are 
incredibly fast by their nature to go low level themselves, so no need to
reinvent wheels here. Nonetheless, the speed-ups you can gain from using Rust can
be incredible. Rust functions should start with *rs_*, and ideally an R wrapper
should exist to use them. Please refer to the (yet to be written) [Why Rust](/docs/why_rust.md)
section.
5. In terms of object-oriented programming, [S7](https://github.com/RConsortium/S7)
provides a way to write very R-like OOP (the methods belong to generics). For 
user-facing key methods and workflows, we recommend using this one, as most R
users will feel very familiar with it and it allows for piping of various 
functions. In  certain cases, you might want to write a more Java/Python-like 
type of classes (i.e., encapsulated OOP), which gives you a lot of control in
terms of public vs. private methods and attributes. We recommend [R6](https://r6.r-lib.org/articles/Introduction.html)
here, but be aware that the average R users might find the R6 classes not very
intuitive. Inheritance can be quite useful in certain cases to abstract out
common generics/methods, but try to avoid deeply layered inheritance where 
possible. This is not the most complex software we are writing here, so there
should be no need for 8+ layers of abstraction.
