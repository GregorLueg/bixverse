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

<img src="/misc/pics/stop_abstracting.png" alt="stop abstracting">
