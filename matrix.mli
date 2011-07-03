(*
Copyright (c) 2011, Svyatoslav Bezgin <bezgin@rcbd.org> All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

 *  Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.
 *  Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.
 *  Neither the name of the Tambov State Technical University nor the names of
    its contributors may be used to endorse or promote products derived from this
    software without specific prior written permission.
            
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*)

val transpose : ('a, 'b, 'c) Batteries.Bigarray.Array2.t -> ('a, 'b, 'c) Batteries.Bigarray.Array2.t

val mx_mult : ('a -> 'b -> 'a) -> ('a -> 'c -> 'b) -> 'a ->
    ('a, 'd, 'e) Batteries.Bigarray.Array2.t -> 
    ('c, 'f, 'g) Batteries.Bigarray.Array2.t -> 
    ('a, 'd, 'e) Batteries.Bigarray.Array2.t

val ( @*. ) :
    (float, 'a, 'b) Batteries.Bigarray.Array2.t ->
    (float, 'c, 'd) Batteries.Bigarray.Array2.t ->
    (float, 'a, 'b) Batteries.Bigarray.Array2.t

val mx_neg : 'a -> ('a -> 'b -> 'b) ->
    ('b, 'c, 'd) Batteries.Bigarray.Array2.t ->
    ('b, 'c, 'd) Batteries.Bigarray.Array2.t

val ( ~@-. ) : (float, 'a, 'b) Batteries.Bigarray.Array2.t -> (float, 'a, 'b) Batteries.Bigarray.Array2.t

val mx_add : ('a -> 'b -> 'a) ->
    ('a, 'c, 'd) Batteries.Bigarray.Array2.t ->
    ('b, 'e, 'f) Batteries.Bigarray.Array2.t ->
    ('a, 'c, 'd) Batteries.Bigarray.Array2.t

val ( @+. ) : (float, 'a, 'b) Batteries.Bigarray.Array2.t ->
    (float, 'c, 'd) Batteries.Bigarray.Array2.t ->
    (float, 'a, 'b) Batteries.Bigarray.Array2.t

val ( @-. ) : (float, 'a, 'b) Batteries.Bigarray.Array2.t ->
    (float, 'c, 'd) Batteries.Bigarray.Array2.t ->
    (float, 'a, 'b) Batteries.Bigarray.Array2.t

val mx_scalar_mul : ('a -> 'b -> 'b) -> 'a ->
    ('b, 'c, 'd) Batteries.Bigarray.Array2.t ->
    ('b, 'c, 'd) Batteries.Bigarray.Array2.t

val ( &*. ) : float -> (float, 'a, 'b) Batteries.Bigarray.Array2.t ->
    (float, 'a, 'b) Batteries.Bigarray.Array2.t

val mx_rk4_step : (float -> (float, 'a, 'b) Batteries.Bigarray.Array2.t -> (float, 'c, 'd) Batteries.Bigarray.Array2.t) ->
    (float, 'a, 'b) Batteries.Bigarray.Array2.t ->
    float -> float -> (float, 'a, 'b) Batteries.Bigarray.Array2.t

val mx_rk4_up : (float -> (float, 'a, 'b) Batteries.Bigarray.Array2.t -> (float, 'c, 'd) Batteries.Bigarray.Array2.t) ->
    (float, 'a, 'b) Batteries.Bigarray.Array2.t ->
    float -> float -> float -> (float * (float, 'a, 'b) Batteries.Bigarray.Array2.t) list

val mx_rk4_down :(float -> (float, 'a, 'b) Batteries.Bigarray.Array2.t -> (float, 'c, 'd) Batteries.Bigarray.Array2.t) ->
    (float, 'a, 'b) Batteries.Bigarray.Array2.t ->
    float -> float -> float -> (float * (float, 'a, 'b) Batteries.Bigarray.Array2.t) list

val mx_dot_product : 'a -> ('a -> 'b -> 'a) -> ('c -> 'd -> 'b) -> 
    ('c, 'e, 'f) Batteries.Bigarray.Array2.t -> 
    ('d, 'g, 'h) Batteries.Bigarray.Array2.t -> 
    'a

val ( ><. ) : (float, 'a, 'b) Batteries.Bigarray.Array2.t ->
    (float, 'c, 'd) Batteries.Bigarray.Array2.t ->
    float

