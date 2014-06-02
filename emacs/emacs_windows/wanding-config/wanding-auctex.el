;; auctex set up
(require 'tex-site)
;; (load "auctex.el" nil t t)
;; (load "preview-latex.el" nil t t)

(add-hook 'LaTeX-mode-hook
	  (lambda ()
	    
	    ;; default to PDF mode
	    ;; pdflatex
	    ;; (TeX-PDF-mode)
	    (setq TeX-PDF-mode t)
	    (outline-minor-mode)
	    (TeX-fold-mode 1); Automatically activate TeX-fold-mode. C-c C-o C-b to fold
	    (add-to-list 'TeX-command-list '("Evince" "evince %s.pdf" TeX-run-silent nil nil))
	    (add-to-list 'TeX-command-list '("View" "evince %s.pdf" TeX-run-silent nil nil))
	    ;; (add-to-list 'TeX-output-view-style
	    ;; 		 (quote ("^pdf$" "." "evince %o %(outpage)")))
	    (setq-default TeX-command-Show "Evince")
	    (setq reftex-toc-split-windows-horizontally t)
	    (local-set-key (kbd "<home>") 'TeX-command-master)
	    ;; (TeX-fold-buffer)
	    )
	  )

(add-hook 'LaTeX-mode-hook 'turn-on-reftex)
