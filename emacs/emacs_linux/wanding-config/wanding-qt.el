;; qt highlight and indentation.

;; (c-add-style "qt-gnu" '("gnu" 
;;                         (c-access-key .
;; 				      "\\<\\(signals\\|public\\|protected\\|private\\|public slots\\|protected slots\\|private slots\\):")
;; 			(c-basic-offset . 4)))


;; (require 'cc-mode)
;; (setq c-C++-access-key "\\<\\(slots\\|signals\\|private\\|protected\\|public\\)\\>[ \t]*[(slots\\|signals)]*[ \t]*:")
;; (font-lock-add-keywords 'c++-mode '(("\\<\\(Q_OBJECT\\|public slots\\|public signals\\|private slots\\|private signals\\|protected slots\\|protected signals\\)\\>" . font-lock-constant-face)))

;; (when (locate-library cc-mode)
;;   (setq c-font-lock-keywords-3
;; 	(append '("signals" "\\(public\\|protected\\|private\\) slots")
;; 		c-font-lock-keywords-3)))
