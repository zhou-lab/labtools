;;Make return automatically indent
;; Lisp (SLIME) interaction

;; (setq inferior-lisp-program "sbcl")
(add-to-list 'load-path "~/.slime")
(setq slime-lisp-implementations
      '((ptools ("/home/wanding/pathway-tools/pathway-tools" "-lisp") :init slime-init-command)))
(require 'slime)
;; (slime-setup '(slime-repl))

;; (setq slime-complete-symbol-function 'slime-fuzzy-complete-symbol)

;; (add-hook 'lisp-mode-hook '(lambda ()
;; 			     (local-set-key (kbd "RET") 'newline-and-indent)))
;; (add-hook 'emacs-lisp-mode-hook (lambda () (local-set-key (kbd "RET") 'newline-and-indent)))