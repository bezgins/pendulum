ocamlfind ocamlc -package batteries,batteries.syntax -syntax camlp4o -thread -linkpkg -c matrix.mli matrix.ml &&
ocamlfind ocamlc -package batteries,batteries.syntax -syntax camlp4o -thread -linkpkg matrix.cmo pendulum.ml -o pendulum
