ocamlfind ocamlopt -package batteries,batteries.syntax -syntax camlp4o -thread -linkpkg -c matrix.mli matrix.ml &&
ocamlfind ocamlopt -package batteries,batteries.syntax -syntax camlp4o -thread -linkpkg matrix.cmx pendulum.ml -o pendulum
