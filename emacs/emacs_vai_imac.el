;; you need to install
;; tabbar
;; yasnippet, ESS, polymode, helm, toc-org, fill-column-indicator
;; other things: ivy+counsel+swiper

(when (>= emacs-major-version 24)
  (require 'package)
  (add-to-list 'package-archives '("melpa" . "http://melpa.milkbox.net/packages/") t)
  ;; (add-to-list 'package-archives '("marmalade" . "http://marmalade-repo.org/packages/"))
	;; (add-to-list 'package-archives '("ELPA", "http://tromey.com/elpa/"))
  (package-initialize)
  )

;;;;;;;;;;;;;;;;;;;;
;;;;;;; fonts ;;;;;;
;;;;;;;;;;;;;;;;;;;;
;; I must say the default emacs font is pretty good to the eyes already
;; (add-to-list 'default-frame-alist '(font . "Inconsolata-12"))
;; (add-to-list 'default-frame-alist '(font . "Source Code Prod 10"))
;; (set-default-font "Source Code Pro 10")
;; (add-to-list 'default-frame-alist '(font . "Source Code Pro-11"))
;; (add-to-list 'default-frame-alist '(font . "Consolas-12"))
;; (add-to-list 'default-frame-alist '(font . "Latin Modern Mono-12"))
;; (add-to-list 'default-frame-alist '(font . "Monaco-12"))
;; (add-to-list 'default-frame-alist '(font . "Iosevka-12:light"))
;; (add-to-list 'default-frame-alist '(font . "Anonymous Pro-13"))
;; (add-to-list 'default-frame-alist '(font . "Anonymous Pro-16"))
(add-to-list 'default-frame-alist '(font . "Menlo-14"))

;; (require 'ess)
(load "ess-site")
auto-mode-alist (append (list '("\\.c$" . c-mode)	
			      '("\\.tex$" . latex-mode)
			      '("\\.S$" . S-mode)
			      '("\\.s$" . S-mode)
			      '("\\.R$" . R-mode)
			      '("\\.r$" . R-mode)
			      '("\\.html$" . html-mode)
            '("\\.emacs" . emacs-lisp-mode)
						)
						auto-mode-alist)

;; for Rmarkdown to find pandoc
;; (setq markdown-command "/Applications/RStudio.app/Contents/MacOS/pandoc")
(setq markdown-command "/usr/local/bin/pandoc")

;; open file in existing frame
(setq ns-pop-up-frames nil)

(setenv "PATH" (concat (getenv "PATH") ":/usr/local/texlive/2014/bin/x86_64-darwin/"))
(setq exec-path (append exec-path '("/usr/local/texlive/2014/bin/x86_64-darwin/")))

;; (load-theme 'ample-zen t)
;; (load-theme 'airline-badwolf t)
;; (load-theme 'tsdh-light t)
;; (load-theme 'leuven t)
;; (load-theme 'zenburn t)
;; (load-theme 'deeper-blue t)
(load-theme 'adwaita t)
;; (load-theme 'dracula t)
;; (setq-default cursor-type 'bar)

(setq-default indent-tabs-mode nil)

(tabbar-mode)

;; (setq yas-snippet-dirs
;;       '~/emacs.d/.elpa/snippets/
;;       )
(yas-global-mode 1)
;; need to link
;; ln -s `rf ~/wzlib/emacs/emacs_linux/wanding-config/snippets` ~/.emacs.d/snippets
;; (yas/load-directory "~/wzlib/emacs/emacs_linux/wanding-config/snippets/")
(add-to-list 'yas-snippet-dirs "~/wzlib/emacs/emacs_linux/wanding-config/snippets")

;; (global-linum-mode t)

(setq tab-width 2)
(setq-default tab-width 2)
(setq sh-basic-offset 2
      sh-indentation 2)

(setq-default indent-tabs-mode nil) ;; indentation cannot insert tabs
(setq org-startup-indented t)
(setq org-startup-folded 'content)
(setq org-src-fontify-natively t)
(setq org-edit-src-content-indentation 0)
(setq org-src-tab-acts-natively t)
(setq org-src-window-setup 'other-window) ;; see C-h v org-src-window-setup, this shows the src code in the other window

;; R support inside org-babel
;; (useful only for export R code-containing org files)
;; (setq org-babel-R-command "/usr/local/bin/R --no-save")

;; TOC in org mode
(if (require 'toc-org nil t)
    (add-hook 'org-mode-hook 'toc-org-mode)
  (warn "toc-org not found"))

(setq mouse-wheel-scroll-amount '(2))
(setq mouse-wheel-progressive-speed nil)
;; (global-visual-line-mode nil)
(set-default 'truncate-lines t)
;; (define-key org-mode-map "\M-q" 'toggle-truncate-lines)
(global-set-key (kbd "M-q") 'toggle-truncate-lines)

(setq word-wrap nil)
(setq mac-command-modifier 'control)

(add-to-list 'default-frame-alist '(height . 50))
(add-to-list 'default-frame-alist '(width . 100))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; fill-column-indicator mode
;; This is not good, see line-move-visual and fci-handle-line-move-visual
;; (define-globalized-minor-mode
;;   global-fci-mode fci-mode (lambda () (fci-mode 1)))
;; (global-fci-mode 0)

(setq-default fill-column 79)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; visualize long lines > 80
;; using whitespace-mode
;; see https://emacs.stackexchange.com/questions/147/how-can-i-get-a-ruler-at-column-80
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; The following solution just highlight lines that exceeds 80
;; (setq-default
;;  whitespace-line-column 80
;;  whitespace-style
;;  '(face lines-tail))
;; (add-hook 'prog-mode-hook #'whitespace-mode)

(setq-default c-default-style "k&r"
              c-basic-offset 3)

(show-paren-mode 1)
(require 'mouse)
(xterm-mouse-mode t)
(defun track-mouse (e)) 
(setq mouse-sel-mode t)
(delete-selection-mode 1)

(global-set-key [mouse-4] (lambda ()
			    (interactive)
			    (scroll-down 1)))
(global-set-key [mouse-5] (lambda ()
			    (interactive)
			    (scroll-up 1)))
(global-set-key (kbd "<f12>") 'save-buffer)
(global-set-key (kbd "<f11>") 'dabbrev-expand)
;; (global-set-key (kbd "<f9>") 'my-switch-to-other-buffer)
(global-set-key (kbd "<f19>") 'my-switch-to-other-buffer)
;; (global-set-key (kbd "<f9>") 'switch-to-buffer)
(global-set-key (kbd "<f18>") 'switch-to-buffer)
(global-set-key (kbd "<f10>") 'yas-reload-all)
(global-set-key (kbd "<f9>") 'yas-describe-tables)

;; copy/paste
(defun zhou-copy-line-or-region ()
  "Copy current line, or region if active."
  (interactive)
  (if (use-region-p)
      (kill-ring-save (region-beginning) (region-end))
    (kill-ring-save (line-beginning-position)
                    (line-end-position)))
  (message "line/region copied"))

(global-set-key (kbd "<f1>") 'zhou-copy-line-or-region) ; copy
(global-set-key (kbd "<f2>") 'yank)           ; paste

;;;;;;;;;;;;;;;;;;;;;;;
;; Wanding functions
;;;;;;;;;;;;;;;;;;;;;;;
(defun my-switch-to-other-buffer ()
  "Switch to other buffer"
  (interactive)
  (switch-to-buffer (other-buffer)))

(defun my-switch-to-previous-buffer ()
  (interactive)
  (switch-to-buffer (other-buffer (current-buffer) 1)))

(defun select-current-line ()
  "Select the current line"
  (interactive)
  (end-of-line) ; move to end of line
  (set-mark (line-beginning-position)))
(global-set-key (kbd "C-'") 'select-current-line)


(defun flush-blank-lines (start end)
  (interactive "r")
  (flush-lines "^\\s-*$" start end nil))

(defun collapse-blank-lines (start end)
  (interactive "r")
  (replace-regexp "^\n\\{2,\\}" "\n" nil start end))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; paragraph-sentence convertion
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; collapse sentences into a paragraph
(defun make-paragraph (start end)
  (interactive "r")
  (replace-regexp "\n\\{2,\\}" " " nil start end))

(global-set-key (kbd "<f6>") 'make-paragraph)

;; split a paragraph into sentences
(defun split-paragraph (start end)
  (interactive "r")
  (replace-regexp "\\. " ".\n\n" nil start end))

(global-set-key (kbd "<f5>") 'split-paragraph)

;; the following is not useful at all, just for learning purpose
(defun remove-vowel ($string &optional $from $to)
  "Remove the following letters: {a e i o u}.

When called interactively, work on current paragraph or text selection.

When called in lisp code, if ξstring is non-nil, returns a changed string.
If ξstring nil, change the text in the region between positions ξfrom ξto."
  (interactive
   (if (use-region-p)
       (list nil (region-beginning) (region-end))
     (let ((bds (bounds-of-thing-at-point 'paragraph)) )
       (list nil (car bds) (cdr bds)) ) ) )

  (let (workOnStringP inputStr outputStr)
    (setq workOnStringP (if $string t nil))
    (setq inputStr (if workOnStringP $string (buffer-substring-no-properties $from $to)))
    (setq outputStr
          (let ((case-fold-search t))
            (replace-regexp-in-string "a\\|e\\|i\\|o\\|u\\|" "" inputStr) )  )

    (if workOnStringP
        outputStr
      (save-excursion
        (delete-region $from $to)
        (goto-char $from)
        (insert outputStr) )) ) )


;; kill whole word implementation
(defcustom brutal-word-regex "[-_A-Za-z0-9]"
  "Regular expression that defines a character
   allowed in a word, i.e., not word boundary. Should
   match only one character at a time."
  :type 'regexp
  :group 'user)

(defun brutally-kill-word ()
  "Kills the whole word (as defined by the
   brutal-word-regex) regardless of where the point is in it."
  (interactive)
  ;; Work only if point is inside of a word
  (if (looking-at brutal-word-regex)
      (save-excursion ; Back up till the word boundary, one char at a time.
	(while (looking-at brutal-word-regex) (backward-char)) 
	;; (message "looking successful")
	;; Search for the whole word ...
	(search-forward-regexp (concat brutal-word-regex "+") nil t)
	;; ... and replace it with empty string
	(replace-match "")
	;; Remove extra spaces if they are there
	(if (eq (char-after) ? ) (just-one-space)))))

(global-set-key (kbd "M-d") 'brutally-kill-word)

(defun rename-file-and-buffer (new-name)
  ;; source: http://steve.yegge.googlepages.com/my-dot-emacs-file
  "Renames both current buffer and file it's visiting to NEW-NAME."
  (interactive "sNew name: ")
  (let ((name (buffer-name))
    (filename (buffer-file-name)))
    (if (not filename)
    (message "Buffer '%s' is not visiting a file!" name)
      (if (get-buffer new-name)
      (message "A buffer named '%s' already exists!" new-name)
    (progn
      (rename-file name new-name 1)
      (rename-buffer new-name)
      (set-visited-file-name new-name)
      (set-buffer-modified-p nil))))))

(defun move-buffer-file (dir)
  "Moves both current buffer and file it's visiting to DIR." 
  (interactive "DNew directory: ")
  (let* ((name (buffer-name))
	 (filename (buffer-file-name))
	 (dir
	  (if (string-match dir "\\(?:/\\|\\\\)$")
	      (substring dir 0 -1) dir))
	 (newname (concat dir "/" name)))
    
    (if (not filename)
	(message "Buffer '%s' is not visiting a file!" name)
      (progn 
	(copy-file filename newname 1) 
	(delete-file filename)
 	(set-visited-file-name newname)
 	(set-buffer-modified-p nil)
 	t)))) 

(defun delete-this-buffer-and-file ()
  "Removes file connected to current buffer and kills buffer.
   with dired-x just do 'C-x C-j D yes RET y RET"
  (interactive)
  (let ((filename (buffer-file-name))
        (buffer (current-buffer))
        (name (buffer-name)))
    (if (not (and filename (file-exists-p filename)))
        (error "Buffer '%s' is not visiting a file!" name)
      (when (yes-or-no-p "Are you sure you want to remove this file? ")
        (delete-file filename)
        (kill-buffer buffer)
        (message "File '%s' successfully removed" filename)))))

(defun jao-copy-line ()
  "Copy current line in the kill ring"
  (interactive)
  (kill-ring-save (line-beginning-position)
		  (line-beginning-position 2))
  (message "line copied"))
(global-set-key (kbd "M-l") 'jao-copy-line)

;; M-x print-to-pdf
(defun print-to-pdf ()
  (interactive)
  (ps-spool-buffer-with-faces)
  (switch-to-buffer "*PostScript*")
  (write-file "/tmp/tmp.ps")
  (kill-buffer "tmp.ps")
  (setq cmd (concat "ps2pdf14 /tmp/tmp.ps " (buffer-name) ".pdf"))
  (shell-command cmd)
  (shell-command "rm /tmp/tmp.ps")
  (message (concat "Saved to:  " (buffer-name) ".pdf"))
  )

(defun move-line-upward ()
  (interactive)
  (kill-region (line-beginning-position) (line-beginning-position 2))
  (goto-char (line-beginning-position 0))
  (yank)
  (goto-char (line-beginning-position 0))
  )


(defun move-line-downward ()
  (interactive)
  (kill-region (line-beginning-position) (line-beginning-position 2))
  (goto-char (line-beginning-position 2))
  (yank)
  (goto-char (line-beginning-position 0))
)


(defun insert-before-line ()
  "insert the content yanked (e.g., copied) before the current line"
  (interactive)
  (let ((pos (point))
        (cur-max (point-max)))
    (beginning-of-line)

    ;; I've changed the order of (yank) and (indent-according-to-mode)
    ;; in order to handle the case when yanked line comes with its own indent
    (yank);; (flush-blank-lines (point))(indent-according-to-mode)
    ;; ;; could be as well changed to simple (newline) it's matter of taste
    ;; ;; and of usage
    ;; (newline-and-indent) 
    ;; back to the original position
    ;; (goto-char (+ pos (- (point-max) cur-max)))
    )
  )
(global-set-key (kbd "C-.") 'insert-before-line)

(defun open-line-above ()
  "Open a line above the line the point is at.
Then move to that line and indent according to mode"
  (interactive)
  (move-beginning-of-line 1)
  (newline)
  (previous-line)
  (indent-according-to-mode))
(global-set-key (kbd "C-o") 'open-line-above)
(global-set-key (kbd "C-<up>") 'open-line-above)

(defun open-line-below ()
  "Open a line below the line the point is at.
Then move to that line and indent according to mode"
  (interactive)
  (move-end-of-line 1)
  (newline)
  (indent-according-to-mode))
(global-set-key (kbd "C-,") 'open-line-below)
;; (global-set-key (kbd "C-m") 'open-line-below)

;; (global-set-key (kbd "C-<down>") 'open-line-below)

(defun reload-dot-emacs()
  "Reload .emacs on the fly"
  (interactive)
  (if(bufferp (get-file-buffer ".emacs"))
      (save-buffer(get-buffer ".emacs")))
  (load-file "~/.emacs")
  (message ".emacs reloaded successfully"))

(defun server-shutdown ()
  "Save buffers, Quit, and Shutdown (kill) server"
  (interactive)
  (save-some-buffers)
  (kill-emacs)
  )

(defun my-put-file-name-on-clipboard ()
  "Put the current file name on the clipboard"
  (interactive)
  (let ((filename (if (equal major-mode 'dired-mode)
                      default-directory
                    (buffer-file-name))))
    (when filename
      (let ((x-select-enable-clipboard t)) (kill-new filename))
      ;; (with-temp-buffer
      ;;   (insert filename)
      ;;   (clipboard-kill-region (point-min) (point-max)))
      (message filename))))

(defun comment-or-uncomment-region-or-line ()
    "Comments or uncomments the region or the current line if there's no active region."
    (interactive)
    (let (beg end)
        (if (region-active-p)
            (setq beg (region-beginning) end (region-end))
            (setq beg (line-beginning-position) end (line-end-position)))
        (comment-or-uncomment-region beg end);))
        (next-line)))
(global-set-key (kbd "C-;") 'comment-or-uncomment-region-or-line)

(defun insert-date ()
  "Insert current date yyyy-mm-dd."
  (interactive)
  (when (region-active-p)
    (delete-region (region-beginning) (region-end) )
    )
  (insert (format-time-string "%Y%m%d_"))
  )

(global-set-key "\C-x\C-d" 'insert-date)

(defun timestamp ()
  (interactive)
  (insert (format-time-string "%Y-%m-%d %H:%M")))

;;;;;;;; UNCLEAR FUNCTIONS ;;;;;;;;;;
(defun insert-and-indent-line-above ()
  "[UNCLEAR FUNC] insert the content before the current line and indent the current line"
  (interactive)
  (push-mark)
  (let* 
    ((ipt (progn (back-to-indentation) (point)))
     (bol (progn (move-beginning-of-line 1) (point)))
     (indent (buffer-substring bol ipt)))
    (newline)
    (previous-line)
    (insert indent)))

;; ;;Place all backup copies of files in a common location
;; (defconst use-backup-dir t)   
;; (setq backup-directory-alist (quote ((".*" . "~/emacs-meta/backups/")))
;;       version-control t                ; Use version numbers for backups
;;       kept-new-versions 16             ; Number of newest versions to keep
;;       kept-old-versions 2              ; Number of oldest versions to keep
;;       delete-old-versions t            ; Ask to delete excess backup versions?
;;       backup-by-copying-when-linked t) ; Copy linked files, don't rename.

(global-set-key (kbd "M-i") 'kill-whole-line)
(global-set-key (kbd "<f16>") 'other-window)
(global-set-key (kbd "<f18>") 'delete-other-windows)
(global-set-key (kbd "<clear>") 'goto-line)
(global-set-key (kbd "M-<up>") 'move-line-upward)
(global-set-key (kbd "M-<down>") 'move-line-downward)
(global-set-key (kbd "<home>") 'beginning-of-line)
(global-set-key (kbd "<end>") 'end-of-line)

;; helm
(setq recentf-max-saved-items 10000)
(global-set-key (kbd "M-r") 'helm-recentf)

(defun insert-file-name ()
  "Insert the full path file name into the current buffer."
  (interactive)
  (insert (buffer-file-name (window-buffer (minibuffer-selected-window)))))
(put 'upcase-region 'disabled nil)

;;;;;;;;;;
;; ESS
;;;;;;;;;;
(require 'ess-r-mode)
;; Note that there are a lot of out-dated information out
;; there on the internet. ESS has changed a lot in terms of
;; the variables for customization. Look at ess-custom.el
;; and C-h v ess-style-alist for the most updated information.
;;
;; To customize, you need to change the default list to OWN
;; all the other style are FIXED!!
;; The default style is RRR. Even changing this thing to
;; DEFAULT won't work. So far only the following works.
(custom-set-variables
 '(ess-own-style-list
   (quote
    ((ess-indent-offset . 4)
     (ess-offset-arguments . prev-line)
     (ess-offset-arguments-newline . prev-line)
     (ess-offset-block . prev-line)
     (ess-offset-continued . straight)
     (ess-align-nested-calls)
     (ess-align-arguments-in-calls)
     (ess-align-continuations-in-calls)
     (ess-align-blocks control-flow)
     (ess-indent-from-lhs)
     (ess-indent-from-chain-start)
     (ess-indent-with-fancy-comments . t)))))

(setq ess-style 'OWN)
(ess-toggle-underscore nil)
(setq inferior-R-program-name "/usr/local/bin/R")

;; (defun rmd-mode ()
;;  "ESS Markdown mode for rmd files"
;;  (interactive)
;;  (require 'poly-R)
;;  (require 'poly-markdown)     
;;  (poly-markdown+r-mode))
(add-to-list 'auto-mode-alist '("\\.Rmd" . poly-markdown+r-mode))
;; the shortcut for R
(defun then_R_operator ()
  "R - %>% operator or 'then' pipe operator"
  (interactive)
  (just-one-space 1)
  (insert "%>%")
  ;; (reindent-then-newline-and-indent)
  (just-one-space 1))
(define-key ess-mode-map (kbd "C-%") 'then_R_operator)
(define-key inferior-ess-mode-map (kbd "C-%") 'then_R_operator)

(setq column-number-mode t)

(defun create-tags (dir-name)
     "Create tags file."
     (interactive "DDirectory: ")
     (eshell-command 
      (format "find %s -type f -name \"*.[ch]\" | etags -" dir-name)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; I don't find ivy mode that compelling compared to helm
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; (ivy-mode 1)
;; (setq ivy-use-virtual-buffers t)
;; (setq enable-recursive-minibuffers t)
;; (global-set-key "\C-s" 'swiper)
;; (global-set-key (kbd "C-c C-r") 'ivy-resume)
;; (global-set-key (kbd "<f6>") 'ivy-resume)
;; (global-set-key (kbd "M-x") 'counsel-M-x)
;; (global-set-key (kbd "C-x C-f") 'counsel-find-file)
;; (global-set-key (kbd "<f1> f") 'counsel-describe-function)
;; (global-set-key (kbd "<f1> v") 'counsel-describe-variable)
;; (global-set-key (kbd "<f1> l") 'counsel-find-library)
;; (global-set-key (kbd "<f2> i") 'counsel-info-lookup-symbol)
;; (global-set-key (kbd "<f2> u") 'counsel-unicode-char)
;; (global-set-key (kbd "C-c g") 'counsel-git)
;; (global-set-key (kbd "C-c j") 'counsel-git-grep)
;; (global-set-key (kbd "C-c k") 'counsel-ag)
;; (global-set-key (kbd "C-x l") 'counsel-locate)
;; (define-key minibuffer-local-map (kbd "C-r") 'counsel-minibuffer-history)
(custom-set-variables
 ;; custom-set-variables was added by Custom.
 ;; If you edit it by hand, you could mess it up, so be careful.
 ;; Your init file should contain only one such instance.
 ;; If there is more than one, they won't work right.
 '(package-selected-packages
   (quote
    (abs-mode markdown-mode yaml-mode yasnippet-snippets yasnippet toc-org tabbar helm fill-column-indicator ess))))
(custom-set-faces
 ;; custom-set-faces was added by Custom.
 ;; If you edit it by hand, you could mess it up, so be careful.
 ;; Your init file should contain only one such instance.
 ;; If there is more than one, they won't work right.
 )
