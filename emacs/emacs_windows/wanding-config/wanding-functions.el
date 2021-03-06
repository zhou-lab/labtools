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
  (insert (format-time-string "%Y_%m_%d_"))
  )

(global-set-key "\C-c\C-d" 'insert-date)

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

