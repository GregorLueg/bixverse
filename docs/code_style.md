# *bixverse package* code style

The general idea is to leverage modern R language features (for example using 
the new [S7 OOP system](https://github.com/RConsortium/S7) to have a more "R 
like" feeling to objects and classes, while avoiding some short falling of 
S3/S4), and use the [rextendr](https://github.com/extendr/rextendr) into Rust 
to make key computations as fast as possible. </br>
Born out of experience, we would like to avoid unnecessary complexity, 
inheritance structures that feel smart but make debugging a huge pain, and the
motto should be at all times **Keep it Simple, Stupid (KISS)**. 
Please, just think about this meme in doubt:

<img src="/misc/pics/stop_abstracting.png" width="350" height="504" alt="stop abstracting">

## *General rules*

1. Most important one: document the code. [roxygen2](https://roxygen2.r-lib.org) 
makes this trivial, so, please, just document your functions. Provide some 
information what parameters mean, what they expect and what is being returned.
R is dynamically typed, so reasoning without going into the source code can be
in times quite painful.
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
data, renaming, plotting go to Rust and use the rextendr interface to make 
computations go *brrrrrr*. Some libraries such as [igraph](https://r.igraph.org)
are incredibly fast by their nature to go low level themselves, so no need to
reinvent wheels here. The speed-ups you can gain from using Rust can be 
incredible. Rust functions should start with *rs_*, and ideally an R wrapper
should exist to use them. Please refer to the 
5. In terms of object-oriented programming, [S7](https://github.com/RConsortium/S7)
provides a way to write very R-like OOP (the methods belong to generics). For 
user-facing key methods and workflows, we recommend using this one, as most R
users will feel very familiar with it and it allows for piping of various 
functions. In  certain cases, you might want to write a more Java/Python-like 
type of classes (i.e., encapsulated OOP), which gives you a lot of control in
terms of public vs. private methods/attributes. We recommend [R6](https://r6.r-lib.org/articles/Introduction.html)
here, but be aware that the average R users might find the R6 classes not very
intuitive.
