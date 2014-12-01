;; for gnus, the following can also be put into .gnus.el
(require 'gnus)
(add-to-list 'gnus-secondary-select-methods '(nnimap "gmail"
                                  (nnimap-address "imap.gmail.com")
                                  (nnimap-server-port 993)
                                  (nnimap-stream ssl)))


(setq user-mail-address "zhouwanding@gmail.com")
(setq user-full-name "Wanding Zhou")
;; (setq gnus-select-method '(nnml ""))
(setq mail-sources '((imap :server "imap.gmail.com"
			   :port 993
			   :user "zhouwanding"
;;			   :password "secret" 
			   :authentication 'login
			   :stream ssl
			   :fetchflag "\\Seen")))



;;(require 'gnus)
;; (add-to-list 'gnus-secondary-select-methods '(nnimap "gmail"
;;                                   (nnimap-address "imap.gmail.com")
;;                                   (nnimap-server-port 993)
;;                                   (nnimap-stream ssl)))

;; (setq gnus-select-method
;;       '(nnimap "gmail"
;; 	       (nnimap-address "imap.gmail.com")
;; 	       (nnimap-server-port 993)
;; 	       (nnimap-stream ssl)))

;; (setq imap-log t)
;; ;; ;; use nnml as backend
;; ;; (setq gnus-select-method '(nnml ""))

;; set nnimap as backend
(setq gnus-select-method '(nnimap "gmail.com"
				  (nnimap-address "imap.gmail.com")
				  (nnimap-stream ssl)))

(setq gnus-read-active-file 'some)

;; Tree view for groups.  I like the organisational feel this has.
(add-hook 'gnus-group-mode-hook 'gnus-topic-mode)

;; Threads!  I hate reading un-threaded email -- especially mailing
;; lists.  This helps a ton!
(setq gnus-summary-thread-gathering-function
      'gnus-gather-threads-by-subject)

;; Also, I prefer to see only the top level message.  If a message has
;; several replies or is part of a thread, only show the first
;; message.  'gnus-thread-ignore-subject' will ignore the subject and
;; look at 'In-Reply-To:' and 'References:' headers.
(setq gnus-thread-hide-subtree t)
(setq gnus-thread-ignore-subject t)

;; the following doesn't work for now
(setq starttls-use-gnutls t)
(setq starttls-gnutls-program "gnutls-cli")
(setq message-send-mail-function 'smtpmail-sent-it
      send-mail-function 'smtpmail-send-it
      smtpmail-starttls-credentials '(("smtp.gmail.com" 587 nil nil))
      smtpmail-auth-credentials '(("smtp.gmail.com" 587 "zhouwanding@gmail.com" nil))
      smtpmail-default-smtp-server "smtp.gmail.com"
      smtpmail-smtp-server "smtp.gmail.com"
      smtpmail-smtp-service 587
      smtpmail-local-domain "rice.edu")
(require 'smtpmail)
