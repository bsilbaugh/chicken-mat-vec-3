;;;
;;; 3 Vectors and 3X3 matrix operations with float elements
;;;
;;; Copyright 2013 Benjamin Silbaugh
;;;
;;; See LICENSE file for redistribution and modification permissions.

(include "mat-vec-3.scm")

(module float (add sub mul div sqroot mx)
  (import scheme chicken)
  (define add fp+)
  (define sub fp-)
  (define mul fp*)
  (define div fp/)
  (define sqroot fpsqrt)
  (define mx fpmax))

(module fp-mat-vec-3 = (mat-vec-3 float))
