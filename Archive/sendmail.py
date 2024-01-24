import argparse
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.application import MIMEApplication


def send_email(subject, to="", attachment=None, body=None):
    # Set up the email server
    smtp_server = 'smtp.gmail.com'
    smtp_port = 587
    sender_email = ''
    sender_password = ''

    # Create a MIMEText object for the email body
    msg = MIMEMultipart()
    msg['From'] = sender_email
    msg['Subject'] = subject

    # Add recipient
    if to:
        msg['To'] = to

    # Add email body
    if body:
        body_text = MIMEText(body, 'plain')
        msg.attach(body_text)

    # Add attachment
    if attachments:
        for attachment in attachments:
            name=attachment.split("\\")[-1]
            with open(attachment, 'rb') as file:
                part = MIMEApplication(file.read())
                part.add_header('Content-Disposition', f'attachment; filename="{name}"')
                msg.attach(part)

    # Send the email
    with smtplib.SMTP(smtp_server, smtp_port) as server:
        server.starttls()
        server.login(sender_email, sender_password)
        server.sendmail(sender_email, to if to else [], msg.as_string())



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Send an email with attachments.")
    parser.add_argument("subject", help="Subject of the email")
    parser.add_argument("-t","--to", help="Recipient's email address")
    parser.add_argument("-a","--attachment", help="Path to the attachment file")
    parser.add_argument("-b","--body", help="Email body")

    args = parser.parse_args()
    attachments = args.attachment.split(",") if args.attachment else []
    # Call the send_email() function with command-line arguments

    if args.to:
        send_email(args.subject, args.to, attachment=attachments, body=args.body)
    else:
        send_email(args.subject, attachment=attachments, body=args.body)




# send_email("testing my app",attachment="D:\Documents\Personal\Projects\BasicFluidSimulation\Final4-1-E1.0-S1.0-C2.5-D1.0-100-5.0-3--0.0001-1000000-100-100-10000-1.0.png")
