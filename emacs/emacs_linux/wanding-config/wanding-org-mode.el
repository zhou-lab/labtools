(defun wanding-org-insert-result ()
  "insert :result: in the beginning of the line"
  (interactive)
  (beginning-of-line)
  (delete-horizontal-space)
  (insert ":result: ")
  (next-line))

(defun wanding-org-insert-data ()
  "insert :data: in the beginning of the line"
  (interactive)
  (beginning-of-line)
  (delete-horizontal-space)
  (insert ":data: ")
  (next-line))

(defun wanding-org-insert-code ()
  "insert :code: in the beginning of the line"
  (interactive)
  (beginning-of-line)
  (insert ":code: ")
  (next-line))

(defun wanding-org-insert-figure ()
  "insert :figure: in the beginning of the line"
  (interactive)
  (beginning-of-line)
  (insert ":figure: ")
  (next-line))

(defun wanding-org-insert-command ()
  "insert :command: in the beginning of the line"
  (interactive)
  (beginning-of-line)
  (insert ":command: ")
  (next-line))

;; (defun zin/org-open-other-frame ()
;;   "Jump to bookmark in another frame. See `bookmark-jump' for more."
;;   (interactive)
;;   (let ((org-link-frame-setup (acons 'file 'find-file-other-frame org-link-frame-setup)))
;;     (org-open-at-point)))
;; (global-set-key (kbd "<mouse-3>") 'zin/org-open-other-frame)

(global-set-key (kbd "M-s r") 'wanding-org-insert-result)
(global-set-key (kbd "M-s d") 'wanding-org-insert-data)
(global-set-key (kbd "M-s c") 'wanding-org-insert-code)
(global-set-key (kbd "M-s f") 'wanding-org-insert-figure)
(global-set-key (kbd "M-s a") 'wanding-org-insert-command)

(global-set-key (kbd "<C-M-return>") 'org-insert-heading-respect-content)

;; (put 'set-goal-column 'disabled nil)
(global-set-key "\C-cl" 'org-store-link)
(global-set-key "\C-ca" 'org-agenda)
(global-set-key (kbd "C-(") 'org-clock-in)
(global-set-key (kbd "C-)") 'org-clock-out)
;; (custom-set-variables
;;   ;; custom-set-variables was added by Custom.
;;   ;; If you edit it by hand, you could mess it up, so be careful.
;;   ;; Your init file should contain only one such instance.
;;   ;; If there is more than one, they won't work right.
;;  '(browse-url-browser-function 'w3m-browse-url)
;;  ;; '(browse-url-browser-function (quote browse-url-generic))
;;  '(browse-url-generic-program "google-chrome")
;;  '(org-agenda-files (quote ("~/Desktop/study/Today.org")))
;;  '(recentf-max-saved-items 100)
;;  '(truncate-partial-width-windows nil))
;; (custom-set-faces
;;   ;; custom-set-faces was added by Custom.
;;   ;; If you edit it by hand, you could mess it up, so be careful.
;;   ;; Your init file should contain only one such instance.
;;   ;; If there is more than one, they won't work right.
;;  )

;; for org mode, use continuation mode of the line instead of trucation mode
;; (add-hook 'org-mode-hook
;; 	  (lambda () (toggle-truncate-lines)))

(add-hook 'org-mode-hook
	  (lambda () (local-set-key (kbd "<C-M-return>") 'org-insert-heading-respect-content)
	    (local-set-key (kbd "C-,") 'open-line-below) ;; override the default org-cycle-agenda
	    ))

(setq org-startup-truncated nil)
;; Org-mode settings
;; (add-to-list 'auto-mode-alist '("\\.org$" . org-mode))
;; (global-font-lock-mode 1)
;; (outline-mode)



;; setup the column that a tag is to be aligned
(setq org-tags-column -50)
(setq org-hide-leading-stars t)
(setq org-agenda-todo-list-sublevels nil)
(setq org-src-fontify-natively t)
;; (setq org-src-fontify-natively nil)
(setq org-todo-keywords
      '((sequence "TODO(t)" "|" "DONE(d)" "CANCELED(c)")))
;; startup visibility
;; (setq org-startup-folded "content")

;; (setq org-pretty-entities t)
(setq org-startup-folded (quote content))
