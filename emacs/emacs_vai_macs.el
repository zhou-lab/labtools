
;; you need to install tabbar, yasnippet, ESS, polymode, helm, toc-org, fill-column-indicator
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
(add-to-list 'default-frame-alist '(font . "Anonymous Pro-13"))

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
;; need to link ~/emacs.d/snippets to ~/wzlib/emacs/emacs_linux/wanding-config/snippets
;; (yas/load-directory "~/wzlib/emacs/emacs_linux/wanding-config/snippets/")
;; (add-to-list 'yas-snippet-dirs "~/wzlib/emacs/emacs_linux/wanding-config/snippets")

;; (global-linum-mode t)

(setq tab-width 2)
(setq-default tab-width 2)
(setq sh-basic-offset 2
      sh-indentation 2)

(setq-default indent-tabs-mode nil)
(setq org-startup-indented t)
(setq org-startup-folded 'content)
(setq org-src-fontify-natively t)
(setq org-edit-src-content-indentation 0)
(setq org-src-tab-acts-natively t)

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
(set-default 'truncate-lines nil)
;; (define-key org-mode-map "\M-q" 'toggle-truncate-lines)
(global-set-key (kbd "M-q") 'toggle-truncate-lines)

(setq word-wrap nil)
(setq mac-command-modifier 'control)

(add-to-list 'default-frame-alist '(height . 50))
(add-to-list 'default-frame-alist '(width . 100))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; fill-column-indicator mode
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define-globalized-minor-mode
  global-fci-mode fci-mode (lambda () (fci-mode 1)))
(global-fci-mode 1)
(set-fill-column 79)
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

;; (setq c-default-style "linux")

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
(global-set-key (kbd "<f13>") 'kill-ring-save)
(global-set-key (kbd "<f14>") 'jao-copy-line)
(global-set-key (kbd "<f16>") 'yank)


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

(custom-set-variables
 ;; custom-set-variables was added by Custom.
 ;; If you edit it by hand, you could mess it up, so be careful.
 ;; Your init file should contain only one such instance.
 ;; If there is more than one, they won't work right.
 '(custom-safe-themes
   (quote
    ("274fa62b00d732d093fc3f120aca1b31a6bb484492f31081c1814a858e25c72e" "c006bc787154c31d5c75e93a54657b4421e0b1a62516644bd25d954239bc9933" "c237299ad2385b4e59767eeb55a0a7e888ebbabec9975be9ae63eec2ff74668d" "f08c2405a7e71e568b784ae0145a86e48e1b4ea8ba33d231a4ad21b52495de5e" "5b86be6c49b591ec04aa238b929f0f4a8d79364653c8482ba69db698645b4de1" "ab8033276aa563bc7373f78aeefef69e1e25083266b44116341f0a8096657463" "40b7687853f3d6921eba3afed50c532bbf4a66959554d32adf1899a675926b2d" "7db66dafe7a65a8a6a403014edb5e53deca2da82279cb8f3f55e4bc336bf48af" "de309af2ced9914b67077eecd0b89412dd9a60c5eb823e5c5ed66170bd4495a7" "de05e8c13f7b8f3f3b9aaee44855c43dcdfd4db8b93c15af2bdcaf2528154ebc" "b69df114abdbbf223e1ad2c98ad1abee04ac2a5070aeb8b7ceefcf00aa5e43f8" "de8fa309eed1effea412533ca5d68ed33770bdf570dcaa458ec21eab219821fd" "3b0f554ddd413e74b82854d78c7c22df6cb4298413f69b514b2884fa84d42f30" "e8a9dfa28c7c3ae126152210e3ccc3707eedae55bdc4b6d3e1bb3a85dfb4e670" "1db337246ebc9c083be0d728f8d20913a0f46edc0a00277746ba411c149d7fe5" "eaf4cb94ad96e1659f9254db8efb799deb1885e97884f8f971ff1e6a4114500a" "356f57a98f35c8ead5a349408cab69f8d4d92baea131e9531611d0d82190fedf" default)))
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
     (ess-indent-with-fancy-comments . t))))
 '(package-selected-packages
   (quote
    (fill-column-indicator dracula-theme toc-org counsel swiper ivy poly-R tabbar-ruler polymode php-mode markdown-mode leuven-theme helm ess-R-object-popup ess auto-yasnippet auctex ample-zen-theme ample-theme airline-themes ahungry-theme))))
(custom-set-faces
 ;; custom-set-faces was added by Custom.
 ;; If you edit it by hand, you could mess it up, so be careful.
 ;; Your init file should contain only one such instance.
 ;; If there is more than one, they won't work right.
 )


(defun insert-file-name ()
  "Insert the full path file name into the current buffer."
  (interactive)
  (insert (buffer-file-name (window-buffer (minibuffer-selected-window)))))
(put 'upcase-region 'disabled nil)

;;;;;;;;;;
;; ESS
;;;;;;;;;;
(setq inferior-R-program-name "/usr/local/bin/R")

(setq ess-style 'OWN)

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
