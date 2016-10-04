## How to input various symbols:

#### ϵ
We use ϵ (lunate epsilon) rather than ε (lower case epsilon). 

This is not the default epsilon for emacs users. You can add something
like this to your .emacs file so that you can easily type an ϵ:
```
(define-key key-translation-map (kbd "<f9> e <right>") (kbd "ϵ"))
```

Then, to insert an ϵ, type `F9 e →`.

#### S''
In LaTeX input mode, type S\prime\prime

