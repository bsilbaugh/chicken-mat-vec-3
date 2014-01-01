;;;
;;; 3 Vector and 3X3 matrix operations with generic arithmetic
;;;
;;; Copyright 2013 Benjamin Silbaugh
;;; 
;;; See LICENSE file for modification and redistribution permissions.

(include "mat-vec-3.scm")

(module num (add sub mul div sqroot mx)
  (import scheme chicken)
  (define add +)
  (define sub -)
  (define mul *)
  (define div /)
  (define sqroot sqrt)
  (define mx max))

(module num-mat-vec-3 = (mat-vec-3 num))
