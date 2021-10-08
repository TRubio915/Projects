using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class IceBall : MonoBehaviour
{
    public float ibSpeed = 5; //Speed of ball
    public Rigidbody2D rb; //RB of ball
    public float fbDamage = 25; //Damage ball does to player
    public Transform hitPoint; //hit point
    public float hitRadius; //Hit radius
    public LayerMask whatIsPlayer; //Player layer mask

    [Header("Touch Hit box")]
    public Vector2 knockBackSpeed;
    public Vector2 knockBackSpeedAttack1;
    public Vector2 knockBackSpeedAttack2;
    public Vector2 knockBackSpeedAttack3;
    public Vector2 touchDamageBotLeft;
    public Vector2 touchDamageTopRight;

    public float
        touchDamageWidth,
        touchDamageHeight;

    public Transform touchDamageCheck;

    private float[] attackDetails;

    public GameObject explodeParticle;


    // Start is called before the first frame update
    void Start()
    {
        attackDetails = new float[2];
        rb.velocity = new Vector2(-1 * ibSpeed, 0); //Move right according to speed
    }


    private void Update()
    {
        touchDamageBotLeft.Set(touchDamageCheck.position.x - touchDamageWidth / 2, touchDamageCheck.position.y - touchDamageHeight / 2);
        touchDamageTopRight.Set(touchDamageCheck.position.x + touchDamageWidth / 2, touchDamageCheck.position.y + touchDamageHeight / 2);

        Collider2D hit = Physics2D.OverlapArea(touchDamageBotLeft, touchDamageTopRight, whatIsPlayer);

        if (hit != null)
        {
            hit.GetComponent<Rigidbody2D>().velocity = Vector2.zero;
            attackDetails[0] = fbDamage / 2;
            attackDetails[1] = transform.position.x;
            hit.SendMessage("Damage", attackDetails);
            Destroy(gameObject);
        }
        Destroy(gameObject, 3);
    }

    public void Damage(float[] attackDetails)
    {
        if(attackDetails[0] > 0)
        {
            DeleteObject();
            Instantiate(explodeParticle, transform.position, Quaternion.identity);
        }
    }


    public void DeleteObject()
    {
        Destroy(gameObject);
    }
    public void OnDrawGizmos()
    {
        Gizmos.color = Color.black;
        Vector2 botLeft = new Vector2(touchDamageCheck.position.x - touchDamageWidth / 2, touchDamageCheck.position.y - touchDamageHeight / 2);
        Vector2 botRight = new Vector2(touchDamageCheck.position.x + touchDamageWidth / 2, touchDamageCheck.position.y - touchDamageHeight / 2);
        Vector2 topLeft = new Vector2(touchDamageCheck.position.x - touchDamageWidth / 2, touchDamageCheck.position.y + touchDamageHeight / 2);
        Vector2 topRight = new Vector2(touchDamageCheck.position.x + touchDamageWidth / 2, touchDamageCheck.position.y + touchDamageHeight / 2);
        Gizmos.DrawLine(botLeft, botRight);
        Gizmos.DrawLine(botRight, topRight);
        Gizmos.DrawLine(topRight, topLeft);
        Gizmos.DrawLine(topLeft, botLeft);
    }

}
